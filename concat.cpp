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


int main(int argc, char **argv)
{
  auto alignments =  std::vector<std::unordered_map<std::string, std::string> > (); 

  if(argc == 1 )
    {
      std::cout << "Usage: " << argv[0] << " file[..]" << std::endl; 
      exit(1); 
    }


  for(int i = 1; i < argc; ++i)
    {
      auto &&in = std::ifstream(std::string(argv[i])); 

      auto aln = std::unordered_map<std::string, std::string>(); 
      std::string line; 
      
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

  // get all the taxa 
  auto taxa = std::unordered_set<std::string>(); 
  for(auto aln : alignments)
    for(auto elem : aln)
      taxa.insert(elem.first); 
  
  // sanity check 
  for(auto aln : alignments)
    {
      auto length = aln.begin()->second.size(); 
      for(auto elem : aln ) 
	assert(length == elem.second.size() ); 
    }

  auto finalAlignment = std::unordered_map<std::string,std::string>(); 
  std::ofstream finalModel("aln.model"); 
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


  // write alignment 
  std::ofstream finalAln("aln.phy"); 
  finalAln << finalAlignment.size() << delim << finalAlignment.begin()->second.size() << std::endl; 
  for(auto &elem : finalAlignment)
    {
      auto name = elem.first; 
      name.erase(
		 std::remove_if(name.begin(), name.end(), 
				[]( int c){return (  c == '>' )  ; }),
		 name.end()); 
      finalAln << name << delim << elem.second << std::endl;       
    }

  return 0; 
}

