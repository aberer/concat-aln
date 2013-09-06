#include <iostream>
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



static void writeAlignment(const std::unordered_map<std::string, std::string>  &finalAlignment,  std::string outfilename, bool isFasta) 
{
  // write alignment 
  std::ofstream finalAln( outfilename); 
  
  if(isFasta)
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
  else 
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
}


int main(int argc, char **argv)
{
  auto alignments =  std::vector<std::unordered_map<std::string, std::string> > (); 

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

  bool isFasta = formatString.compare("FASTA") == 0 ; 
  if(not isFasta)
    {
      bool isPhylip = formatString.compare("PHYLIP") == 0; 
      if(not isPhylip)
	{
	  std::cout << "format string (first argument) must be either PHYLIP or FASTA" << std::endl; 
	  exit(-1); 
	}
    }


  
  for(int i = 3; i < argc; ++i)
    {
      auto &&in = std::ifstream(std::string(argv[i])); 

      auto aln = std::unordered_map<std::string, std::string>(); 
      std::string line; 
      
      std::cout << "parsing alignment "  << argv[i] << std::endl;
      
      auto lastElem = aln.begin(); 
      while(std::getline(in, line))
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
      
      alignments.push_back(aln);
    }

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

  auto finalAlignment = std::unordered_map<std::string,std::string>(); 
  std::ofstream finalModel(id + ".model"); 
  nat pos = 1; 
  nat ctr = 1; 
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

      // write model file 
      nat newPos = pos + aln.begin()->second.size(); 
      auto partName = basename(argv[ctr])      ; 
      finalModel << "DNA, " <<  partName << "=" << pos << "-" << newPos -1  << std::endl ; 
      ++ctr; 
      pos = newPos ;
    }


  writeAlignment(finalAlignment, id + ".phy", isFasta); 

  return 0; 
}

