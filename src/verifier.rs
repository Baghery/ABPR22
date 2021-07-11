use ark_ec::msm::FixedBaseMSM;
use ark_ec::{AffineCurve, PairingEngine, ProjectiveCurve};
use ark_ff::{PrimeField, to_bytes, Field,};

use super::{PreparedVerifyingKey, Proof, VerifyingKey};

use ark_relations::r1cs::{Result as R1CSResult, SynthesisError};

use core::ops::{AddAssign, Neg};
use blake2::{Blake2b, Digest};
use ark_std::{vec::Vec};

/// Prepare the verifying key `vk` for use in proof verification.
pub fn prepare_verifying_key<E: PairingEngine>(vk: &VerifyingKey<E>) -> PreparedVerifyingKey<E> {
    PreparedVerifyingKey {
        vk: vk.clone(),
        gamma_g2_neg_pc: vk.gamma_g2.neg().into(),
    }
}

/// Prepare proof inputs for use with [`verify_proof_with_prepared_inputs`], wrt the prepared
/// verification key `pvk` and instance public inputs.
pub fn prepare_inputs<E: PairingEngine>(
    pvk: &PreparedVerifyingKey<E>,
    public_inputs: &[E::Fr],
) -> R1CSResult<E::G1Projective> {
    if (public_inputs.len() + 1) != pvk.vk.gamma_abc_g1.len() {
        return Err(SynthesisError::MalformedVerifyingKey);
    }

    let mut g_ic = pvk.vk.gamma_abc_g1[0].into_projective();
    for (i, b) in public_inputs.iter().zip(pvk.vk.gamma_abc_g1.iter().skip(1)) {
        g_ic.add_assign(&b.mul(i.into_repr()));
    }

    Ok(g_ic)
}

/// Verify a proof `proof` against the prepared verification key `pvk` and prepared public
/// inputs. This should be preferred over [`verify_proof`] if the instance's public inputs are
/// known in advance.
pub fn verify_proof_with_prepared_inputs<E: PairingEngine>(
    pvk: &PreparedVerifyingKey<E>,
    proof: &Proof<E>,
    prepared_inputs: &E::G1Projective,
) -> R1CSResult<bool> {
    let qap = E::miller_loop(
        [
            (proof.a.into(), proof.b.into()),
            (
                prepared_inputs.into_affine().into(),
                pvk.gamma_g2_neg_pc.clone(),
            ),
            (proof.c.into(), proof.delta_prime.clone().neg().into()),
        ]
        .iter(),
    );

    let test = E::final_exponentiation(&qap).ok_or(SynthesisError::UnexpectedIdentity)?;

    let hash = Blake2b::new()
    .chain(to_bytes!(&proof.a).unwrap())
    .chain(to_bytes!(&proof.b).unwrap())
    .chain(to_bytes!(&proof.c).unwrap())
    .chain(to_bytes!(&proof.delta_prime).unwrap());
    let mut output = [0u8; 64];
    output.copy_from_slice(&hash.finalize());

    let m_fr = E::Fr::from_le_bytes_mod_order(&output);

    let tmp1 = pvk.vk.delta_g2.clone().mul(m_fr);
    let delta_prime_delta_m = proof.delta_prime + tmp1.into_affine();
    let test2 = E::pairing(proof.d, delta_prime_delta_m);
    
    Ok((test == pvk.vk.alpha_g1_beta_g2) && (test2 == pvk.vk.zt_gt.pow(m_fr.into()) * pvk.vk.kappa_zt_gt))
}

/// Verify a proof `proof` against the prepared verification key `pvk`,
/// with respect to the instance `public_inputs`.
pub fn verify_proof<E: PairingEngine>(
    pvk: &PreparedVerifyingKey<E>,
    proof: &Proof<E>,
    public_inputs: &[E::Fr],
) -> R1CSResult<bool> {
    let prepared_inputs = prepare_inputs(pvk, public_inputs)?;
    verify_proof_with_prepared_inputs(pvk, proof, &prepared_inputs)
}






/// Verify a vector of proofs `proofs` against the prepared verification key `pvk` and prepared public
/// inputs. This should be preferred over [`verify_proof`] if the instance's public inputs are
/// known in advance.
pub fn vec_verify_proof_with_prepared_inputs<E: PairingEngine>(
    pvk: &PreparedVerifyingKey<E>,
    proofs: &Vec<Proof<E>>,
    prepared_inputs: &Vec<E::G1Projective>,
) -> R1CSResult<bool> {
    let mut m_fr: Vec<E::Fr> = Vec::new();
    let mut size_proofs = 0usize;
    for proof in proofs.iter() {
        let hash = Blake2b::new()
        .chain(to_bytes!(&proof.a).unwrap())
        .chain(to_bytes!(&proof.b).unwrap())
        .chain(to_bytes!(&proof.c).unwrap())
        .chain(to_bytes!(&proof.delta_prime).unwrap());
        let mut output = [0u8; 64];
        output.copy_from_slice(&hash.finalize());
        m_fr.push(E::Fr::from_le_bytes_mod_order(&output));
        size_proofs +=1 ;
    }

    
    let scalar_bits = E::Fr::size_in_bits();
    let delta_g2_window = FixedBaseMSM::get_mul_window_size(size_proofs);
    let delta_g2_table =
        FixedBaseMSM::get_window_table::<E::G2Projective>(scalar_bits, delta_g2_window, pvk.vk.delta_g2.into_projective());
    let elem_g2 =
        FixedBaseMSM::multi_scalar_mul::<E::G2Projective>(scalar_bits, delta_g2_window, &delta_g2_table, &m_fr);


    let num_powers = scalar_bits;
    let mut powers_of_2 = Vec::with_capacity(num_powers as usize);

    let mut p = pvk.vk.zt_gt;
    powers_of_2.push(p);
    for _ in 1..num_powers {
        p.square_in_place();
        powers_of_2.push(p);
    }
    
   
    let result = 
    elem_g2.iter().zip(proofs).zip(m_fr).zip(prepared_inputs).map(|(((x, y), z),w)|  
    (E::pairing(y.d,y.delta_prime + x.into_affine()) == 
                        E::Fqk::pow_with_table(&powers_of_2[..], z.into()).unwrap() * pvk.vk.kappa_zt_gt)
    &&
    (E::final_exponentiation(&E::miller_loop(
        [
            (y.a.into(), y.b.into()),
            (
                w.into_affine().into(),
                pvk.gamma_g2_neg_pc.clone(),
            ),
            (y.c.into(), y.delta_prime.neg().into()),
        ]
        .iter(),
    )).
    unwrap() == pvk.vk.alpha_g1_beta_g2)).
    fold(true, |total, next| {total && next});
    
     
    //println!("result is {:?}", result);
    
    Ok(result)
}

/// Verify a vector of proof `proofs` against the prepared verification key `pvk`,
/// with respect to the instance `public_inputs`.
pub fn vec_verify_proof<E: PairingEngine>(
    vk: &VerifyingKey<E>,
    proofs: &Vec<Proof<E>>,
    public_inputs: &Vec<Vec<E::Fr>>,
) -> R1CSResult<bool> {
    let mut prepared_inputs: Vec<_> = Vec::new();
    for (_,pub_input) in public_inputs.iter().enumerate(){
        let pvk = prepare_verifying_key(vk);
        prepared_inputs.push(prepare_inputs(&pvk, pub_input)?);
    }
    let pvk = prepare_verifying_key(vk);
    vec_verify_proof_with_prepared_inputs(&pvk, proofs, &prepared_inputs)
}