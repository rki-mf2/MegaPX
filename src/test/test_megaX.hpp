#ifndef BUILD_FM_HPP_
#define BUILD_FM_HPP_


#include <filesystem>

#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/fm_index/bi_fm_index.hpp>
#include <seqan3/search/search.hpp>
 
#include <fstream> // for writing/reading files
 
#include <cereal/archives/binary.hpp> // for storing/loading indices via cereal

using namespace seqan3::literals;
using namespace std::string_literals;


void buildIndex(const std::string& inputFastaFile){

    std::vector<seqan3::dna4> text{"ACGTAGC"_dna4};
 
    auto hashes = text | seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{3}});
    seqan3::debug_stream << hashes << '\n'; // [6,27,44,50,9]
 
    seqan3::debug_stream << (text | seqan3::views::kmer_hash(seqan3::ungapped{3})) << '\n'; // [6,27,44,50,9]
 
    seqan3::debug_stream << (text | seqan3::views::kmer_hash(0b101_shape)) << '\n'; // [2,7,8,14,1]

    seqan3::debug_stream << typeid(hashes).name() << std::endl;


    using traits_type = seqan3::sequence_file_input_default_traits_aa;
    seqan3::sequence_file_input<traits_type> fin{inputFastaFile};

    std::vector<seqan3::aa27_vector> sequences; 

    for (auto & record : fin)
    {
        //seqan3::debug_stream << "ID:  " << record.id() << '\n';
        //seqan3::debug_stream << "SEQ: " << record.sequence() << '\n';
        sequences.push_back(record.sequence());
        // a quality field also exists, but is not printed, because we know it's empty for FASTA files.
    }

    std::cout << sequences.size() << std::endl;

    std::cout << "Start building normal index...."<< std::endl;
    seqan3::fm_index sequences_index_normal{sequences};
    std::cout << "Finished building normal index!...."<< std::endl;

    {
        std::ofstream os{"normal_index.idx", std::ios::binary};
        cereal::BinaryOutputArchive oarchive{os};
        oarchive(sequences_index_normal);
    }

    std::cout << "Start building bidirectional index...."<< std::endl;
    seqan3::bi_fm_index sequences_index_bi{sequences};
    std::cout << "Finished building bidirectional index!...."<< std::endl;

    {
        std::ofstream os{"bi_index.idx", std::ios::binary};
        cereal::BinaryOutputArchive oarchive{os};
        oarchive(sequences_index_normal);
    }

    

}
#endif /* BUILD_FM_HPP_ */