
use std::time::Instant;
use colored::*;

mod error;
mod fasta;
mod matrices;
mod options;
mod result;
mod spmeasure;
mod weighting;

fn main()
{
	println!( "\nCalculating conservation score a site in a MSA using Sum-of-pairs measure.\n" );
	println!( "Valdar, William SJ. \"Scoring residue conservation.\" Proteins: structure, function, and bioinformatics 48.2 (2002)." );

	/* Elapsed time : Start */
	let start = Instant::now();

	/* Set options. */
	let opts = options::Options::new();
	opts.show_parameter();

	/* Read an input file and get FASTA information. */
	let mut data = fasta::Fasta::new();
	data.read_fasta_info( &( opts.input ) );

	/* Check whether the input file is correct FASTA format. */
	data.check_fasta_info( &( opts.tolerate ) );

	/* Get site information as Vec<String>. */
	data.get_site_list();

	/*
	println!( "\nInputfile content :\n" );
	for i in 0 .. ( data.seq_list ).len() {
		println!( "Title    {} : {}", i + 1, ( data.title_list )[ i ] );
		println!( "Sequence {} : {}", i + 1, ( data.seq_list )[ i ] );
	}
	*/

	/*
	println!( "\nSite content :\n" );
	for i in 0 .. ( data.site_list ).len() {
		println!( "Site {} : {}", i + 1, ( data.site_list )[ i ] );
	}
	*/

	/* Sequence weighting. */
	let weight_list : Vec<f64> = weighting::seq_weight
	(
		&( data.seq_list  ),
		&( data.site_list ),
		&( opts.weight    )
	);

	/*
	println!( "\nSequence weighting :\n" );
	for i in 0 .. weight_list.len() {
		println!( "Weighing factor of Sequence {} : {}", i + 1, weight_list[ i ] );
	}
	*/

	let cons_valdar01_list : Vec<f64> = spmeasure::sum_of_pairs
	(
		&( data.site_list ),
		&weight_list,
		&( opts.matrix    )
	);

	/*
	for i in 0 .. cons_valdar01_list.len() {
		println!( "Sum-of-pairs measure site {} : {:.3}", i + 1, cons_valdar01_list[ i ] );
	}
	*/

	/* Show result */
	result::show_result
	(
		&( data.site_list ),
		&cons_valdar01_list,
		&( opts.colorize  )
	);

	/* Save result */
	result::save_result
	(
		&( data.site_list ),
		&cons_valdar01_list,
		&( opts.output    )
	);

	println!( "{}", "\nProgram completed !!!\n".green() );

	/* Elapsed time : End */
	let end = start.elapsed();
	println!( "Total elapsed time : {:?}", end );
}
