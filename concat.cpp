#include <iostream>
#include <cstring>
#include <fstream>
#include <libgen.h>
#include <unordered_set>
#include <cassert>
#include <algorithm>
#include <vector>
#include <unordered_map>
#include <string>
#include <sstream>

typedef unsigned int nat; 
const char delim = '\t'; 

enum class FORMAT
{
  FASTA = 1, 
    PHYLIP = 2
}; 


static void writeModelFile(const std::vector<std::unordered_map<std::string, std::string>> &alignments, std::string outfile, const std::vector<std::string > &names)
{
  assert(names.size() == alignments.size() ); 

  std::ofstream finalModel(outfile); 

  nat pos = 1; 
  nat ctr = 0; 
  for(auto aln : alignments)
    {
      nat newPos = pos + aln.begin()->second.size(); 
      auto strCpy = strdup(names.at(ctr).c_str()); 
      auto partName = std::string(basename(strCpy) ); 
      finalModel << "DNA, " <<  partName << "=" << pos << "-" << newPos -1  << std::endl ; 
      ++ctr; 
      pos = newPos ;
      free(strCpy); 
    }
}


// just a quick check 
static FORMAT guessInputFormat(std::string fileName) 
{
  std::ifstream fh(fileName); 
  auto line = std::string{}; 
  std::getline(fh, line); 

  bool found = false; 
  found |= line.find('>') != std::string::npos; 
  
  std::getline(fh,line); 
  found |= line.find('>') != std::string::npos; 
  

  if(not found )
    return FORMAT::PHYLIP; 
  else 
    return FORMAT::FASTA;   
}




auto parseAlignments(const std::vector<std::string> &fileNames)
-> std::vector<std::unordered_map<std::string, std::string> >
{
  auto result = std::vector<std::unordered_map<std::string, std::string> >{}; 

  for(auto fileName : fileNames)
    {
      auto &&in = std::ifstream(fileName); 
      
      auto format = guessInputFormat(fileName);

      nat numTax = 0; 
      nat seqLen = 0; 
      if(format == FORMAT::PHYLIP) 
	{
	  auto line = std::string{}; 
	  std::getline(in, line); 
	  
	  auto &&iss = std::istringstream(line); 
	  iss >> numTax; 
	  iss >> seqLen; 	  
	}

      auto aln = std::unordered_map<std::string, std::string>(); 
      std::string line; 
      
      std::cout << "parsing alignment "  << fileName << std::endl;
      
      auto lastElem = aln.begin(); 
      while(std::getline(in, line))
	{
	  switch(format)
	    {
	    case FORMAT::FASTA: 
	      {	  
		line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end());

		if(line[0]== '>')
		  {
		    assert(aln.find(line) == aln.end()); 
		    aln.emplace(line, ""); 
		    lastElem = aln.find(line); 
		  }
		else 
		  {
		    lastElem->second += line; 
		  }
	      }
	      break; 
	    case FORMAT::PHYLIP:
	      {
		// critical NOTICE: i am only parsing very simplisticly composed input files   
		auto &&iss = std::istringstream(line); 
		auto taxon = std::string{}; 
		auto seq = std::string{}; 
		iss >> taxon; 
		taxon = ">" + taxon; 
		std::getline(iss, seq); 
		seq.erase(std::remove_if(seq.begin(), seq.end(), isspace), seq.end());
		assert(aln.find(taxon) == aln.end()); 
		assert(seq.size( )== seqLen); 
		aln[taxon] =  seq; 
	      }
	      break; 
	    default : 
	      assert(0); 
	    }
	}      

      if(format == FORMAT::PHYLIP )
	{
	  assert(aln.size() == numTax); 
	}
      
      result.push_back(aln);
    }

  return result; 
}



static std::unordered_map<std::string,std::string> concatenateAlignment(const std::unordered_set<std::string> &taxa,const std::vector<std::unordered_map<std::string, std::string> > &alignments )
{
  auto finalAlignment =   std::unordered_map<std::string,std::string> (); 
  for(auto aln : alignments)
    {
      auto length = aln.begin()->second.size(); 
      for(auto taxon : taxa)
	{
	  auto found = aln.find(taxon); 
	  if( found == aln.end() )
	    finalAlignment[taxon] += std::string(length,'-'); 
	  else 	    
	    finalAlignment[taxon] += found->second; 
	}
    }
  return finalAlignment; 
}


static void writeAlignment(const std::unordered_map<std::string, std::string>  &finalAlignment,  std::string outfilename, FORMAT format) 
{
  // write alignment 
  std::ofstream finalAln( outfilename); 

  switch(format)
    {
    case FORMAT::FASTA: 
      {
	for(auto &elem : finalAlignment)
	  {
	    auto name = elem.first; 
	    name.erase( std::remove_if(name.begin(), name.end(), 
				       []( int c){return (  c == '>' )  ; }),
			name.end()); 
	    finalAln << ">" << name << "\n" ; 
	    finalAln << elem.second << std::endl; 
	  }
      }
      break; 
    case FORMAT::PHYLIP: 
      {
	finalAln << finalAlignment.size() << delim << finalAlignment.begin()->second.size() << std::endl; 
	for(auto &elem : finalAlignment)
	  {
	    auto name = elem.first; 
	    name.erase( std::remove_if(name.begin(), name.end(), 
				       []( int c){return (  c == '>' )  ; }),
			name.end()); 
	    finalAln << name << delim << elem.second << std::endl;       
	  }
      }
      break; 
    default : assert(0); 
    }
}


int main(int argc, char **argv)
{

  if(argc < 3 )
    {
      std::cout << "Usage: " << argv[0] << "\tFASTA|PHYLIP id file[..]"
		<< "\n\twhere FASTA|PHYLIP indicates the format of the output file\n"
		<< "\tand id is an id that is used for creating the output files"
		<< std::endl; 
      exit(1); 
    }

  auto formatString = std::string(argv[1]); 
  auto id = std::string(argv[2]); 
  
  auto outputFormat = FORMAT::FASTA; 
  bool isFasta = formatString.compare("FASTA") == 0 ; 
  if(not isFasta)
    {
      bool isPhylip = formatString.compare("PHYLIP") == 0; 
      if(not isPhylip)
	{
	  std::cout << "format string (first argument) must be either PHYLIP or FASTA" << std::endl; 
	  exit(-1); 
	}
      outputFormat = FORMAT::PHYLIP; 
    }

  // generate names 
  auto fileNames = std::vector<std::string>{}; 
  for(int i =  3 ; i < argc; ++i)
    fileNames.push_back(std::string(argv[i]));


  auto alignments = parseAlignments(fileNames); 
  std::cout << "parsed alignments" << std::endl;

  // get all the taxa 
  auto taxa = std::unordered_set<std::string>(); 
  for(auto aln : alignments)
    for(auto elem : aln)
      taxa.insert(elem.first); 
  std::cout << "extracted taxa " << std::endl; 

  // sanity check 
  for(auto aln : alignments)
    {
      auto length = aln.begin()->second.size(); 
      for(auto elem : aln ) 
	assert(length == elem.second.size() ); 
    }
  std::cout << "checked" << std::endl;


  auto finalAlignment = concatenateAlignment(taxa, alignments); 

  writeModelFile(alignments, id + ".model", fileNames); 
  writeAlignment(finalAlignment, id + ".phy", outputFormat); 

  return 0; 
}

