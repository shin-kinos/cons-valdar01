
use std::collections::HashMap;
use std::f64;

use crate::error;
use crate::matrices;

/* Substitution scoring matrix */
static mut M : Vec<f64> = Vec::new();

pub fn sum_of_pairs( site_list : &Vec<String>, weight_list : &Vec<f64>, arg_m : &String ) -> Vec<f64>
{
	let num_site : usize = ( *site_list ).len();

	/* 
	 * Amino acid index for picking the substitution scoring matrix elements.
	 * aa_index = Amino acid index
	 * aa_list  = Amino acid order of substitution scoring matrices
	 */
	let mut aa_index : HashMap<char, usize> = HashMap::new();
	let     aa_list  : Vec<char> = "ARNDCQEGHILKMFPSTWYV".chars().collect();
	for i in 0 .. 20 {
		aa_index.insert( aa_list[ i ], i );
	}
	//println!( "aa_index : {:?}", aa_index );

	/* Make a substitution scoring matrix. */
	unsafe { 
		/* Define a  substitution scoring matrix. */
		M = matrices::define_matrix( arg_m );
		M.shrink_to_fit();
		//println!( "{:?}", M );
		//println!( "\nMaxima of the matrix : {:.3}",   M.iter().fold( 0.0 / 0.0, | m, v | v.max( m ) ) );
		//println!(   "Minima of the matrix : {:.3}\n", M.iter().fold( 0.0 / 0.0, | m, v | v.min( m ) ) );

		/* Check the amino acid index in the substitution scoring matrix. */
		/*
		for i in aa_list.iter() {
			for j in aa_list.iter() {
				println!( "M[ {}, {} ] = {}", *i, *j, M[ ( aa_index )[ i ] * 20 + ( aa_index )[ j ] ] );
			}
		}
		*/

		/* Normalise the scoring matrix based on Karlin-like method */
		M = normalize_matrix( /* &aa_list, &aa_index */ );
		/*
		for a in aa_list.iter() {
			for b in aa_list.iter() {
				println!( "Normalized M[ {}, {} ] = {:.3}", *a, *b, M[ aa_index[ a ] * 20 + aa_index[ b ] ] );
			}
		}
		*/
		//println!( "\nMaxima of the matrix : {:.3}",   M.iter().fold( 0.0 / 0.0, | m, v | v.max( m ) ) );
		//println!(   "Minima of the matrix : {:.3}\n", M.iter().fold( 0.0 / 0.0, | m, v | v.min( m ) ) );

		/* Check whether the matrix is diagonal. */
		check_matrix_diag();

	}

	/* 
	 * Calculate λ ( a value for bounding the conservation score from 0 to 1 ). 
	 * λ = 1 / ( Σ Wi * Wj, i < j )
	 */
	let lambda : f64 = calc_lambda( &weight_list );
	println!( "lambda = {:.3}", lambda );

	/* Calculate the conservation score Valdar01. */
	let mut sp_score_list : Vec<f64> = vec![ 0.0; num_site ];
	for i in 0 .. num_site {
		let mut sp_score : f64 = calc_sp( &( *site_list )[ i ], &weight_list, &aa_index );
		sp_score = lambda * sp_score; 
		//println!( "sp_score of site {} : {} ", i + 1, sp_score );
		sp_score_list[ i ] = sp_score;
	}

	sp_score_list
}

fn normalize_matrix( /* aa_list : &Vec<char>, aa_index : &HashMap<char, usize> */ ) -> Vec<f64>
{
	let mut norm_matrix : Vec<f64> = vec![ 0.0; 400 ];

	/*
	 * Normalise the scoring matrix based on Linear method.
	 * It ensures that the elements of the matrix are bounded from 0 to 1.
	 * mmax = Maxima of the matrix
	 * mmin = Minima of the matrix
	 */
	unsafe {
		let mmax : f64 = M.iter().fold( 0.0 / 0.0, | m, v | v.max( m ) );
		let mmin : f64 = M.iter().fold( 0.0 / 0.0, | m, v | v.min( m ) );
		//println!( "{}, {}", mmax, mmin );
		for a in 0 .. 20 {
			for b in 0 .. 20 {
				//let mat_aa : f64 = M[ a * 20 + a ];
				//let mat_bb : f64 = M[ b * 20 + b ];
				let mat_ab : f64 = M[ a * 20 + b ];
				//norm_matrix[ a * 20 + b ] = mat_ab / ( mat_aa * mat_bb ).sqrt();
				norm_matrix[ a * 20 + b ] = ( mat_ab - mmin ) / ( mmax - mmin );
			}
		}
	}

	/* Check the size of the matrix. */
	if norm_matrix.len() != 400 {
		error::error_bomb( "mat_not_20*20" );
	}

	norm_matrix
}

fn calc_lambda( weight_list : &Vec<f64>) -> f64
{
	let mut lambda : f64 = 0.0;
	let weight_list_len : usize = ( *weight_list ).len();

	for i in 0 .. ( weight_list_len -1 ) {
		for j in ( i + 1 ) .. weight_list_len {
			lambda += ( *weight_list )[ i ] * ( *weight_list )[ j ];
		}
	}

	lambda = 1.0 / lambda;

	lambda
}

fn calc_sp( site : &String, weight_list : &Vec<f64>, aa_index : &HashMap<char, usize> ) -> f64
{
	let     char_list : Vec<char> = ( *site ).chars().collect();
	let mut sp_score  : f64       = 0.0;
	let     site_len  : usize     = char_list.len();

	unsafe {
		/*
		 * Calculate residue conservation using Sum-of-pairs measure.
		 * If ( i = gap ) or ( j = gap ), 0 is given as gap penalty.
		 * i             = Site i
		 * j             = Site j
		 * site_len      = Length of a site
		 * char_list     = Vec<char> of the site (&String)
		 * aa_index[]    = The order of amino acids
		 * weight_list[] = weighting factors
		 * M[]           = Normalized scoring matrix
		 * sp_score      = Conservation score
		 */
		for i in 0 .. ( site_len - 1 ) {
			for j in ( i + 1 ) .. site_len {
				if char_list[ i ] == '-' {
					//println!( "M[ -, {} ] = {:.3}", char_list[ j ], 0.0 );
					sp_score += 0.0;
				} else if char_list[ j ] == '-' {
					//println!( "M[ {}, - ] = {:.3}", char_list[ i ], 0.0 );
					sp_score += 0.0;
				} else {
					let mat_order : usize = ( *aa_index )[ &char_list[ i ] ] * 20 + ( *aa_index )[ &char_list[ j ] ];
					//println!( "M[ {}, {} ] = {:.3}", char_list[ i ], char_list[ j ], ( *weight_list )[ i ] * ( *weight_list )[ j ] * M[ mat_order ] );
					sp_score += ( *weight_list )[ i ] * ( *weight_list )[ j ] * M[ mat_order ];
				}
			}
		}
	}

	//sp_score *= 2.0 / ( site_len * ( site_len - 1 ) ) as f64;

	sp_score
}

fn check_matrix_diag()
{
	unsafe {
		for a in 0 .. 20 {
			for b in 0 .. 20 {
				if M[ a * 20 + b ] != M[ b * 20 + a ] {
					error::error_bomb( "mat_not_diag" );
				}
			}
		}

	}

}
