#include <iostream>
#include <fstream>
#include <map>

#include <seqan/bam_io.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/align.h>
#include <seqan/graph_msa.h> 

#include <basics.h>
#include <variant.h>
#include <merge_split_alignments.h>
#include <evaluate_final_mapping.h>
      
using namespace std;
using namespace seqan;
 
int const BUFFER = 400; // size of flanking region for breakpoints

int const CIGAR_MAX_LEN = 250;

int main(int argc, char ** argv)
{
    if (argc != 3)
    {
        std::cerr << "USAGE: final.sam original.vcf";
        return 1;
    }
 
    string vcf_filename{argv[2]};

    // dorukb
    string bam_filename{argv[1]};

    ifstream vcf_file_in(vcf_filename.c_str());
    ofstream vcf_file_out((vcf_filename + ".polished.vcf").c_str());
    BamFileIn bamfileIn;



    if (!vcf_file_in.is_open())
    {
        std::cerr << "ERROR: Could not open " << argv[2] << std::endl;
        return 1;
    }

    if (!open(bamfileIn, argv[1]))
    { 
        std::cerr << "ERROR: Could not open " << argv[1] << std::endl;
        return 1;
    }


    // dorukb
    BamFileOut bamfileOut(context(bamfileIn), (bam_filename + ".merged_alignments.sam").c_str()); 
    int records_omitted_forbeinglong = 0;

    // -------------------------------------------------------------------------
    // Read in Variants
    // -------------------------------------------------------------------------
    map<string, Variant> variant_map;
    std::string line{"#"};

    while(line[0] == '#')                    // skip header
        getline(vcf_file_in, line);
 
    Variant tmp_var(line);
    variant_map.insert(make_pair(tmp_var.id, tmp_var)); // insert first variant

    while (getline(vcf_file_in, line))          // get the rest of the variants
    {
        tmp_var = Variant(line);
        variant_map.insert(make_pair(tmp_var.id, tmp_var));
    }
  
    // -------------------------------------------------------------------------
    // Merge supplementary alignments before evaluating
    // -------------------------------------------------------------------------
  
    BamHeader header;
    readHeader(header, bamfileIn);

    //dorukb
    writeHeader(bamfileOut, header);
    // empty file must be check here otherwise the first read record will fail
    if (atEnd(bamfileIn))
        return 0; 
 
    BamAlignmentRecord record;
    vector<BamAlignmentRecord> record_group; // will contain all records with the same read name
 
    // BAM file must be sorted by name
    while (!atEnd(bamfileIn))
    {
        readRecord(record, bamfileIn);

        if (!record_group.empty() &&
            (record_group[record_group.size() - 1]).qName != record.qName)
        {
            BamAlignmentRecord merged_record = merge_record_group(record_group);

            // Evaluate alignment and ouput new polished variant
            // -----------------------------------------------------------------
            Variant polished_variant = evaluate_alignment(merged_record, variant_map);

            // dorukb 
            bool has_hardclip = false;
            for (auto ce : merged_record.cigar)
            {
                if (ce.operation == 'H')
                {
                    has_hardclip = true;
                    break;
                }
                
            }
                             
            merged_record.qual = "*";
            if (!has_hardclip)
            {
                writeRecord(bamfileOut, merged_record);
            }
            else
            {
                records_omitted_forbeinglong += 1;
            }
            
            if (polished_variant.sv_type != SV_TYPE::UNKOWN)
            {
                polished_variant.write(vcf_file_out);
            }

            record_group.clear();
            record_group.push_back(record);
        }
        else
        {
            record_group.push_back(record);
        }
    }
    //process last group
    BamAlignmentRecord merged_record = merge_record_group(record_group);
    

    // dorukb
    //std::cout << merged_record.qual  << ": is original" << std::endl;
    merged_record.qual = "*"; 
    bool has_hardclip = false;
    unsigned int fragment_length = compute_fragment_length(merged_record.cigar);
    for (auto ce : merged_record.cigar){
        if (ce.operation == 'H')
        {
            has_hardclip = true;
            break;
        }
    }
                     
    if (!has_hardclip)
    {
        writeRecord(bamfileOut, merged_record);
    }
    else
    {
        records_omitted_forbeinglong += 1;
    }

 
 
    Variant polished_variant = evaluate_alignment(merged_record, variant_map);
    if (polished_variant.sv_type != SV_TYPE::UNKOWN)
    {
        polished_variant.write(vcf_file_out);
    }

    cout << "Done evaluating. Results are in file "
         << vcf_filename + ".polished.vcf" << endl;

    cout << "Records ommited due to having too long of a CIGAR string: "
        << records_omitted_forbeinglong << endl;

    return 0;
}
