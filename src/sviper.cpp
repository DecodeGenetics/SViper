#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>

#include <seqan/arg_parse.h>
#include <seqan/bam_io.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/align.h>
#include <seqan/graph_msa.h>

#include <basics.h>
#include <config.h>
#include <merge_split_alignments.h>
#include <polishing.h>
#include <variant.h>

using namespace std;
using namespace seqan;
  
struct CmdOptions 
{
    bool verbose{false};
    bool veryVerbose{false};
    int flanking_region{400}; // size of flanking region for breakpoints
    double ont_coverage{40}; // ont coverage 
    string long_read_file_name;
    string short_read_file_name;
    string candidate_file_name;
    string out_fa_file_name;
    string reference_file_name;
    string log_file_name;
};  
 
ArgumentParser::ParseResult
parseCommandLine(CmdOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("SViper");
    setVersion(parser, "2.0.0");

    addOption(parser, seqan::ArgParseOption(
        "c", "candidate-vcf",
        "A structural variant vcf file (with e.g. <DEL> tags), containing the potential variant sites to be looked at.",
        seqan::ArgParseArgument::INPUT_FILE, "VCF_FILE"));

    addOption(parser, seqan::ArgParseOption(
        "s", "short-read-bam",
        "The indexed bam file containing short used for polishing at variant sites.",
        seqan::ArgParseArgument::INPUT_FILE, "BAM_FILE"));

    addOption(parser, seqan::ArgParseOption(
        "l", "long-read-bam",
        "The indexed bam file containing long reads to be polished at variant sites.",
        seqan::ArgParseArgument::INPUT_FILE, "BAM_FILE"));

    addOption(parser, seqan::ArgParseOption(
        "r", "reference",
        "The indexed (fai) reference file.",
        seqan::ArgParseArgument::INPUT_FILE, "FA_FILE"));

    addOption(parser, seqan::ArgParseOption(
        "k", "flanking-region",
        "The flanking region in bp's around a breakpoint to be considered for polishing",
        seqan::ArgParseArgument::INTEGER, "INT"));
    //added by dorukb 
    addOption(parser, seqan::ArgParseOption(
        "x", "ont-coverage",
        "The coverage of the ONT bam file",
        seqan::ArgParseArgument::DOUBLE, "DOUBLE"));

    addOption(parser, seqan::ArgParseOption(
        "o", "output-fa",
        "A filename for the output file. NOTE: The current output is a fasta file, that contains the "
        "polished sequences for each variant. Since the final realignment is not part of this tool yet,"
        "The user must map the fasta file with a mapper of his choice (e.g. minimap2) and then call"
        " evaluate_final_alignment.",
        seqan::ArgParseArgument::INPUT_FILE, "FA_FILE"));

    addOption(parser, seqan::ArgParseOption(
        "g", "log-file",
        "A filename for the log file.",
        seqan::ArgParseArgument::INPUT_FILE, "TXT_FILE"));

    addOption(parser, seqan::ArgParseOption(
        "v", "verbose",
        "Turn on detailed information about the process."));

    // addOption(parser, seqan::ArgParseOption(
    //     "vv", "very-verbose",
    //     "Turn on detailed information about the process and write out intermediate results."));

    setRequired(parser, "c");
    setRequired(parser, "l");
    setRequired(parser, "s");
    setRequired(parser, "r");

    setMinValue(parser, "k", "50");
    setMaxValue(parser, "k", "1000");
    setDefaultValue(parser, "k", "400");


    // added by dorukb
    // TODO: 
    setDefaultValue(parser, "x", "40");

    setDefaultValue(parser, "g", "polishing.log");

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    // Extract option values.
    getOptionValue(options.candidate_file_name, parser, "candidate-vcf");
    getOptionValue(options.long_read_file_name, parser, "long-read-bam");
    getOptionValue(options.short_read_file_name, parser, "short-read-bam");
    getOptionValue(options.reference_file_name, parser, "reference");
    getOptionValue(options.out_fa_file_name, parser, "output-fa");
    getOptionValue(options.log_file_name, parser, "log-file");
    getOptionValue(options.flanking_region, parser, "flanking-region");
    getOptionValue(options.ont_coverage, parser, "ont-coverage");
    options.verbose = isSet(parser, "verbose");
    // options.veryVerbose = isSet(parser, "very-verbose");

    //if (options.veryVerbose)
    //    options.verbose = true;

    if (options.out_fa_file_name.empty())
        options.out_fa_file_name = options.candidate_file_name + "polished.vcf";

    return seqan::ArgumentParser::PARSE_OK;
}

int main(int argc, char const ** argv)
{
    // Parse Command Line Arguments
    // -------------------------------------------------------------------------
    CmdOptions options;
    ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    // Check files
    // -------------------------------------------------------------------------
    ifstream input_vcf;           // The candidate variants to polish
    // ofstream output_vcf;          // The polished variant as output
    FaiIndex faiIndex;            // The reference sequence fasta index
    BamHeader long_read_header;   // The bam header object needed to fill bam context
    BamHeader short_read_header;  // The bam header object needed to fill bam context
    BamFileIn long_read_bam;      // The long read bam file
    BamFileIn short_read_bam;     // THe short read bam file
    BamIndex<Bai> long_read_bai;  // The bam index to the long read bam file
    BamIndex<Bai> short_read_bai; // The bam index to the short read bam file
    SeqFileOut final_fa;          // The current output fasta file to be mapped

    if (!open_file_success(input_vcf, options.candidate_file_name.c_str()) ||
        /*!open_file_success(output_vcf, options.out_fa_file_name.c_str()) || */ // no direct vcf output available yet
        !open_file_success(long_read_bam, options.long_read_file_name.c_str()) ||
        !open_file_success(short_read_bam, options.short_read_file_name.c_str()) ||
        !open_file_success(long_read_bai, (options.long_read_file_name + ".bai").c_str()) ||
        !open_file_success(short_read_bai, (options.short_read_file_name + ".bai").c_str()) ||
        !open_file_success(faiIndex, options.reference_file_name.c_str()) ||
        !open_file_success(final_fa, (options.out_fa_file_name + ".fa").c_str()) ||
        !open_file_success(log_file, options.log_file_name.c_str()))
        return 1;

    std::cout << "======================================================================" << std::endl
              << "START polishing variants in of file " << options.candidate_file_name << std::endl
              << "======================================================================" << std::endl;
    log_file  << "======================================================================" << std::endl
              << "START polishing variants in of file " << options.candidate_file_name << std::endl
              << "======================================================================" << std::endl;

    // Read bam file header
    // -------------------------------------------------------------------------
    // This is neccessary for the bam file context (containing reference ids and
    // more) is correctly initialized.
    readHeader(long_read_header, long_read_bam);
    readHeader(short_read_header, short_read_bam);

    // Scip VCF header
    // -------------------------------------------------------------------------
    std::string line{"#"}; // TODO:: use seqan vcf parser instead (needs to be extended)
    while (line[0] == '#' && getline(input_vcf, line)) {}

    // Polish variants
    // -------------------------------------------------------------------------
    do // per Variant
    {
        auto start = std::chrono::high_resolution_clock::now(); // capture time per variant
        Variant var{line};

        if (var.alt_seq != "<DEL>" && var.alt_seq != "<INS>")
        {
            std::cout << "[ SKIP ] Variant " << var.id << " at " << var.ref_chrom << ":"
                      << var.ref_pos << " of type " << var.alt_seq
                      << " because this type is currently not supported." << std::endl;
            log_file  << "----------------------------------------------------------------------" << std::endl
                      << " SKIP Variant " << var.id << " at " << var.ref_chrom << ":" << var.ref_pos << " " << var.alt_seq << " L:" << var.sv_length << std::endl
                      << "----------------------------------------------------------------------" << std::endl;
            continue;
        }
        if (var.sv_length > 1000000)
        {
            std::cout << "[ SKIP ] Variant " << var.id << " at " << var.ref_chrom << ":"
                      << var.ref_pos << " of type " << var.alt_seq
                      << " of length " << var.sv_length << " because it is too long." << std::endl;
            log_file  << "----------------------------------------------------------------------" << std::endl
                      << " SKIP too long Variant " << var.ref_chrom << ":" << var.ref_pos << " " << var.alt_seq << " L:" << var.sv_length << std::endl
                      << "----------------------------------------------------------------------" << std::endl;
            continue;
        }
  
        log_file << "----------------------------------------------------------------------" << std::endl
                 << " PROCESS Variant " << var.id << " at " << var.ref_chrom << ":" << var.ref_pos << " " << var.alt_seq << " L:" << var.sv_length << std::endl
                 << "----------------------------------------------------------------------" << std::endl;
  
        // Extract long reads
        // ---------------------------------------------------------------------

        // flag to tell that the view_bam had to extract too many reads there for illumina.
        bool too_many_illumina_reads = false;


        std::vector<seqan::BamAlignmentRecord> ont_reads;
 
        // extract overlapping the start breakpoint +-50 bp's
        view_bam(ont_reads, long_read_bam, long_read_bai, var.ref_chrom, max(0, var.ref_pos - 50), var.ref_pos + 50, true, too_many_illumina_reads);
        
        std::cout << "ONT read size after start BP only: " <<  ont_reads.size() << " at ref pos: " 
        << var.ref_chrom << ":"  << var.ref_pos << std::endl;
 
        // extract overlapping the end breakpoint +-50 bp's
        view_bam(ont_reads, long_read_bam, long_read_bai, var.ref_chrom, max(0, var.ref_pos_end - 50), var.ref_pos_end + 50, true, too_many_illumina_reads);

        std::cout << "ONT read size after both BPs: " <<  ont_reads.size() << " at ref pos: "  
        << var.ref_chrom << ":" << var.ref_pos << std::endl;


        //store the unique names positions of the reads.
        std::set< std::tuple<std::string, unsigned int, bool > > ont_read_traits;
        for (auto & o_read : ont_reads)
        {
            
            bool non_primary_status = (hasFlagSupplementary(o_read) || hasFlagSecondary(o_read));
            ont_read_traits.insert(std::make_tuple(seqan::toCString(o_read.qName), o_read.beginPos, non_primary_status));
            //std::cout << seqan::toCString(o_read.qName) << " " << o_read.beginPos << " " << non_primary_status << std::endl;
        }


        std::map<std::string, int> ont_read_name_counts;
        for (auto & r_traits : ont_read_traits)
        {
            auto key = std::get<0>(r_traits);
            auto it = ont_read_name_counts.find(key);
            if (it == ont_read_name_counts.end())
            {
                ont_read_name_counts[key] = 1;
            }
            else
            {
                ont_read_name_counts[key] = ont_read_name_counts.at(key) + 1;
            }
            //std::cout << seqan::toCString(o_read.qName) << std::endl;
        }


        int max_read_occurence = 0;
        std::string max_occurence_read = "";

        int too_many_same_read_counts = 3;
        int total_unique_ont_reads = 0;


        for (auto & element : ont_read_name_counts)
        { 
            total_unique_ont_reads += 1;
            auto val = element.second;
            if (val >= max_read_occurence) 
            {   max_read_occurence = val;
                max_occurence_read = element.first;
                // std::cout << element.first << std::endl;;
                // std::cout << val << std::endl; 
            }
        }

        std::cout <<  "max occurence is: " << max_read_occurence << std::endl;
        std::cout << "read name is: " << max_occurence_read << std::endl;
        log_file  <<  "max occurence is: " << max_read_occurence << std::endl;
        log_file  << "read name is: " << max_occurence_read << std::endl;
        if (max_read_occurence >= too_many_same_read_counts)
        {
            std::cout << "marked for removal." << std::endl;
            log_file  << "[ ERRORTOOREPETITIVE ] "
                      << "Current ONT coverage is: " << ont_reads.size() << ". "
                      << var.ref_chrom << ":" << max(0, var.ref_pos - 50) << "-"
                      << var.ref_pos + 50 << " or "
                      << var.ref_chrom << ":" << max(0, var.ref_pos_end - 50) << "-"
                      << var.ref_pos_end + 50 << "or " << std::endl;
            std::cout << "[ ERRORTOOREPETITIVE ] "
                      << "Current ONT coverage is: " << ont_reads.size() << ". "
                      << var.ref_chrom << ":" << max(0, var.ref_pos - 50) << "-"
                      << var.ref_pos + 50 << " or "
                      << var.ref_chrom << ":" << max(0, var.ref_pos_end - 50) << "-"
                      << var.ref_pos_end + 50 << "or " << std::endl;
            continue;
        }



        std::cout << "read names size: " << total_unique_ont_reads << std::endl;

        if (ont_reads.size() == 0)
        {
            log_file  << "[ ERROR1 ] No long reads in reference region "
                      << var.ref_chrom << ":" << max(0, var.ref_pos - 50) << "-"
                      << var.ref_pos + 50 << " or "
                      << var.ref_chrom << ":" << max(0, var.ref_pos_end - 50) << "-"
                      << var.ref_pos_end + 50 << "or " << std::endl;
            std::cout << "[ ERROR1 ] No long reads in reference region "
                      << var.ref_chrom << ":" << max(0, var.ref_pos - 50) << "-"
                      << var.ref_pos + 50 << " or "
                      << var.ref_chrom << ":" << max(0, var.ref_pos_end - 50) << "-"
                      << var.ref_pos_end + 50 << "or " << std::endl;
            continue;
        }
        else if (total_unique_ont_reads > (3*options.ont_coverage) )
        {
            log_file  << "[ ERRORX ] Too many long reads (> 3 times the coverage = "<< (3*options.ont_coverage) <<" ) in the region."
                      << "Current ONT coverage is: " << ont_reads.size() << ". "
                      << var.ref_chrom << ":" << max(0, var.ref_pos - 50) << "-"
                      << var.ref_pos + 50 << " or "
                      << var.ref_chrom << ":" << max(0, var.ref_pos_end - 50) << "-"
                      << var.ref_pos_end + 50 << "or " << std::endl;
            std::cout <<  "[ ERRORX ] Too many long reads (> 3 times the coverage = "<< (3*options.ont_coverage) <<" ) in the region."
                      << "Current ONT coverage is: " << ont_reads.size() << ". "
                      << var.ref_chrom << ":" << max(0, var.ref_pos - 50) << "-"
                      << var.ref_pos + 50 << " or "
                      << var.ref_chrom << ":" << max(0, var.ref_pos_end - 50) << "-"
                      << var.ref_pos_end + 50 << "or " << std::endl;
            continue;
        }

        log_file << "--- Extracted " << ont_reads.size() << " long read(s)" << std::endl;

        // Merge supplementary alignments to primary
        // ---------------------------------------------------------------------
        
        // Hard filter on N number of reads used. Let N be 100.
        // std::sort(ont_reads.begin(), ont_reads.end(), bamRecordMapQGreater());
        // std::cout << "sorting is done: " << ont_reads.size() << std::endl;
        
        // int max_ont_read_tobeused = 100;
        // int total_ont_used = ont_reads.size();
        // int tmp_ont_size = std::min(total_ont_used, max_ont_read_tobeused);
        

        // std::vector<seqan::BamAlignmentRecord> tmp_ont_reads;
        // for (int i = 0; i < tmp_ont_size; i++)
        // {
        //     tmp_ont_reads.push_back(ont_reads[i]);   
        // }
        // ont_reads.clear();
        // for (int i = 0; i < tmp_ont_reads.size(); i++)
        // {
        //     ont_reads.push_back(tmp_ont_reads[i]);
        // }


        // ont_reads.resize(std::min((int)ont_reads.size(), max_ont_read_tobeused));
        //std::cout << "subsetting is done " << ont_reads.size() << std::endl;
        //std::cout << "len of ont_reads: " << ont_reads.size() << std::endl;
        std::cout << "total ont_reads: " << ont_reads.size() << std::endl;

        std::sort(ont_reads.begin(), ont_reads.end(), bamRecordNameLess());
        std::cout << "sort successful. " << std::endl;
        ont_reads = merge_alignments(ont_reads); // next to merging this will also get rid of duplicated reads
       std::cout << "merge successful. " << std::endl;
        log_file << "--- After merging " << ont_reads.size() << " read(s) remain(s)." << std::endl;

        // Search for supporting reads
        // ---------------------------------------------------------------------
        vector<BamAlignmentRecord> supporting_records;
        // TODO check if var is not empty!

        log_file << "--- Searchin in (reference) region ["
                 << (int)(var.ref_pos - DEV_POS * var.sv_length) << "-"
                 << (int)(var.ref_pos + var.sv_length + DEV_POS * var.sv_length) << "]"
                 << " for a variant of type " << var.alt_seq
                 << " of length " << (int)(var.sv_length - DEV_SIZE * var.sv_length) << "-"
                 << (int)(var.sv_length + DEV_SIZE * var.sv_length) << " bp's" << std::endl;

        for (auto const & rec : ont_reads)
            if (record_supports_variant(rec, var))
                supporting_records.push_back(rec);

        if (supporting_records.size() == 0)
        {
            std::cout << "[ ERROR2 ] No supporting reads for a " << var.alt_seq
                      << " in region " << var.ref_chrom << ":"
                      << max(0, var.ref_pos - 50) << "-" << var.ref_pos_end + 50
                      << std::endl;
            continue;
        }

        log_file << "--- After searching for variant " << supporting_records.size()
                 << " supporting read(s) remain." << std::endl;

        // Compute reference length, start and end position of region of interest
        // ---------------------------------------------------------------------
        unsigned ref_fai_idx = 0;
        if (!seqan::getIdByName(ref_fai_idx, faiIndex, var.ref_chrom))
        {
            log_file << "[ ERROR ]: FAI index has no entry for reference name"
                     << var.ref_chrom << std::endl;
            continue;
        }

        int ref_length = seqan::sequenceLength(faiIndex, ref_fai_idx);
        int ref_region_start = std::max(0, var.ref_pos - options.flanking_region);
        int ref_region_end   = std::min(ref_length, var.ref_pos_end + options.flanking_region);

        log_file << "Reference region " << var.ref_chrom << ":" << ref_region_start << "-" << ref_region_end << std::endl;

        SEQAN_ASSERT_LEQ(ref_region_start, ref_region_end);

        // Crop fasta sequence of each supporting read for consensus
        // ---------------------------------------------------------------------
        log_file << "--- Cropping long reads with a buffer of +-" << options.flanking_region << " around variants." << endl;

        StringSet<Dna5String> supporting_sequences;
        std::vector<seqan::BamAlignmentRecord>::size_type maximum_long_reads = 5;

        // sort records such that the highest quality ones are chosen first
        std::sort(supporting_records.begin(), supporting_records.end(), bamRecordMapQGreater());

        #ifndef NDEBUG
                std::cout << "supporting_records.size(): " << supporting_records.size() << std::endl;
        #endif


        bool region_is_incorrect = false;
        bool region_is_toolarge = false;
        for (unsigned i = 0; i < std::min(maximum_long_reads, supporting_records.size()); ++i)
        {
            auto region = get_read_region_boundaries(supporting_records[i], ref_region_start, ref_region_end);

            int32_t expected_length{10*options.flanking_region};
            if (var.sv_type == SV_TYPE::INS)
                expected_length += (2*var.sv_length);

            // added by dorukb as a temporary fix.
            if (get<0>(region) > get<1>(region))
            {
                #ifndef NDEBUG
                std::cout << "get<0>(region): " << get<0>(region) << std::endl;
                std::cout << "get<1>(region): " << get<1>(region) << std::endl;
                #endif
                region_is_incorrect = true;
                break;
            }

            Dna5String reg = seqan::infix(supporting_records[i].seq, get<0>(region), get<1>(region));

            // added by dorukb
            // arbitrary. Filters too long cropped regions. Potentially erroneous and time consuming.
            //else if ((get<1>(region) - get<0>(region)) > (options.flanking_region *10))
            if (static_cast<int32_t>(length(reg)) > expected_length)
            {
                #ifndef NDEBUG
                std::cout << "get<0>(region): " << get<0>(region) << std::endl;
                std::cout << "get<1>(region): " << get<1>(region) << std::endl;
                #endif
                log_file << "------ Region: [" << get<0>(region) << "-" << get<1>(region)
                     << "] Length:" << length(reg) << " and len diff: " 
                     << abs(static_cast<int32_t>(length(reg)) - expected_length) 
                     << " Qual:" << supporting_records[i].mapQ
                     << " Name: "<< supporting_records[i].qName 
                     << " Expected length: " <<  expected_length << endl;
                region_is_toolarge = true;
                //break;
                 ++maximum_long_reads;
                continue; // do not use under or oversized region
            }


            appendValue(supporting_sequences, reg);

            log_file << "------ Region: [" << get<0>(region) << "-" << get<1>(region)
                     << "] Length:" << length(reg) << " and len diff: " 
                     << abs(static_cast<int32_t>(length(reg)) - expected_length) 
                     << " Qual:" << supporting_records[i].mapQ
                     << " Name: "<< supporting_records[i].qName << endl;
        }

        //added by dorukb. As a temporary fix to the incorrect region returned. Debug it.
        if (region_is_incorrect)
        {
            log_file  << "[ ERROR_REGION ] get<0>(region) > get<1>(region) at "
                      << var.ref_chrom << ":" << var.ref_pos << std::endl;
            // std::cout << "[ ERROR_REGION ] get<0>(region) > get<1>(region) at "
            //            << var.ref_chrom << ":" << var.ref_pos << std::endl;
            //var.filter = "FAIL3";
            //#pragma omp critical
            //log_file << localLog.str() << std::endl;
            continue;
        }
        //added by dorukb. As a temporary fix to the incorrect region returned. Debug it.
        // if (region_is_toolarge)
        // {
        //     log_file  << "[ ERROR_REGION ] cropped region in the ONT is too large at "
        //               << var.ref_chrom << ":" << var.ref_pos << std::endl;
        //     std::cout << "[ ERROR_REGION ] cropped region in the ONT is too large at "
        //               << var.ref_chrom << ":" << var.ref_pos   << std::endl;
        //     continue;  
        // }

        if (length(supporting_sequences) == 0)
        {
            log_file << "ERROR3: No fitting regions for a " << var.alt_seq
                     << " in region " << var.ref_chrom <<  ":" << var.ref_pos << std::endl;
            //var.filter = "FAIL3";
            //#pragma omp critical
            //log_file << localLog.str() << std::endl;
            continue;
        }


        // Build consensus of supporting read regions
        // ---------------------------------------------------------------------
        vector<double> mapping_qualities;
        mapping_qualities.resize(supporting_records.size());
        for (unsigned i = 0; i < supporting_records.size(); ++i)
            mapping_qualities[i] = (supporting_records[i]).mapQ;

        Dna5String cns = build_consensus(supporting_sequences, mapping_qualities);

        log_file << "--- Built a consensus with a MSA of length " << length(cns) << "." << endl;

        // Extract short reads in region
        // ---------------------------------------------------------------------
        StringSet<Dna5QString> short_reads_1; // reads (first in pair)
        StringSet<Dna5QString> short_reads_2; // mates (second in pair)

        vector<BamAlignmentRecord> short_reads;
        // If the breakpoints are farther apart then illumina-read-length + 2 * flanking-reagion,
        // then ectract reads for each break point seperately.
        if (ref_region_end - ref_region_start > options.flanking_region * 2 + 150) // 150 = illumina length
        {
            // extract reads left of the start of the variant [start-flanking_region, start+flanking_region]
            unsigned e = min(ref_length, var.ref_pos + options.flanking_region);
            view_bam(short_reads, short_read_bam, short_read_bai, var.ref_chrom, ref_region_start, e, false, too_many_illumina_reads);
            // and right of the end of the variant [end-flanking_region, end+flanking_region]
            unsigned s = max(0, var.ref_pos_end - options.flanking_region);
            view_bam(short_reads, short_read_bam, short_read_bai, var.ref_chrom, s, ref_region_end, false, too_many_illumina_reads);
        }
        else 
        {
            // extract reads left of the start of the variant [start-flanking_region, start]
            view_bam(short_reads, short_read_bam, short_read_bai, var.ref_chrom, ref_region_start, ref_region_end, false, too_many_illumina_reads);
        }
        std::cout << "total short_reads: " << short_reads.size() << std::endl;
 
        if (too_many_illumina_reads == true)
        {
             log_file << "ERRORN: too many illumina reads "
                     << " in region " << var.ref_chrom <<  ":" << var.ref_pos << std::endl;
            //var.filter = "FAIL3";
            //#pragma omp critical
            //log_file << localLog.str() << std::endl;
            continue;
        }


        
        if (short_reads.size() < 20)
        {
            std::cout << "[ ERROR3 ] Not enough short reads (only " << short_reads.size()
                      << ") for variant of type " << var.alt_seq
                      << " in region " << var.ref_chrom << ":" << ref_region_start
                      << "-" << ref_region_end << std::endl;
            continue;
        }

        records_to_read_pairs(short_reads_1, short_reads_2, short_reads, short_read_bam, short_read_bai);

        if (length(short_reads_1) > 250) // ~400 reads -> 60x coverage is enough
            subsample(short_reads_1, short_reads_2, 250); // sub sample short reads to 250 pairs

        log_file << "--- Extracted " << short_reads.size() << " short reads that come in "
                 << length(short_reads_1) << " pairs (proper or dummy pairs)." << std::endl;

        // Flank consensus sequence
        // ---------------------------------------------------------------------
        // Before polishing, append a reference flank of length 150 (illumina
        // read length) to the conesnsus such that the reads find a high quality
        // anchor for mapping.
        String<char> id = "consensus";
        Dna5String flanked_consensus = append_ref_flanks(cns, faiIndex,
                                                         ref_fai_idx, ref_length,
                                                         ref_region_start, ref_region_end, 150);

        // Polish flanked consensus sequence with short reads
        // ---------------------------------------------------------------------
        SViperConfig config{}; // default
        config.verbose = options.verbose;
        compute_baseQ_stats(config, short_reads_1, short_reads_2); //TODO:: return wualities and assign to cinfig outside

        log_file << "--- Short read base qualities: avg=" << config.baseQ_mean
                 << " stdev=" << config.baseQ_std << "." << std::endl;

        Dna5String polished_ref = polish_to_perfection(short_reads_1,
                                                       short_reads_2,
                                                       flanked_consensus, config);

        log_file << "DONE POLISHING: Total of "
                 << config.substituted_bases << " substituted, "
                 << config.deleted_bases     << " deleted and " 
                 << config.inserted_bases    << " inserted bases. "
                 << config.rounds            << " rounds."
                 << std::endl << std::endl;
 
        // Flank polished sequence 
        // ---------------------------------------------------------------------
        // Append large reference flank for better mapping results of long reads.
        Dna5String final_sequence = append_ref_flanks(polished_ref,
                                                      faiIndex, ref_fai_idx,
                                                      ref_length,
                                                      ref_region_start - 150,
                                                      ref_region_end + 150,
                                                      4850);

        // Print polished and flanked sequence to file
        // ---------------------------------------------------------------------
        std::string read_identifier = (string("final") +
                                       ":" + var.ref_chrom +
                                       ":" + to_string(var.ref_pos) +
                                       ":" + var.id); // this name is important for evaluating

        writeRecord(final_fa, read_identifier, final_sequence);

        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "[ SUCCESS ] Polished variant " << var.id << " at " << var.ref_chrom << ":"
                  << var.ref_pos << "\t" << var.alt_seq << "\t[ "
                  << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << " s]"
                  << std::endl;
    }
    while (getline(input_vcf, line));

    std::cout << "======================================================================" << std::endl
              << "                                 DONE"  << std::endl
              << "======================================================================" << std::endl;
    log_file  << "======================================================================" << std::endl
              << "                                 DONE"  << std::endl
              << "======================================================================" << std::endl;

    return 0;
}
