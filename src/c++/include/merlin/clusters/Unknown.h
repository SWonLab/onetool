////////////////////////////////////////////////////////////////////// 
// clusters/Unknown.h 
// (c) 2000-2007 Goncalo Abecasis
// 
// This file is distributed as part of the MERLIN source code package   
// and may not be redistributed in any form, without prior written    
// permission from the author. Permission is granted for you to       
// modify this file for your own personal use, but modified versions  
// must retain this copyright notice and must not be distributed.     
// 
// Permission is granted for you to use this file to compile MERLIN.    
// 
// All computer programs have bugs. Use this file at your own risk.   
// 
// Tuesday December 18, 2007
// 
 
#ifndef __UNKNOWN_H__
#define __UNKNOWN_H__

#include "merlin/libsrc/IntArray.h"
#include "merlin/libsrc/StringMap.h"

#define U_GRAPH      0
#define U_MISSING    1
#define U_INVALID    2

class TUnknown
   {
   public:
      int type;
      int marker;

	  TUnknown()
         {
         marker = -1;
         type = U_INVALID;
         }

	  virtual ~TUnknown() {}
   };

class MissingAllele : public TUnknown
   {
   public:
      int haplotype;

      MissingAllele()
         { type = U_MISSING; };

      virtual ~MissingAllele() {}
   };

class AlleleGraph : public TUnknown
   {
   public:
      IntArray haplotypes;
      IntArray alleles;

      AlleleGraph()
         { type = U_GRAPH; }

      virtual ~AlleleGraph() {}

      void Append(int haplo, int allele1, int allele2)
         {
         haplotypes.Push(haplo);
         alleles.Push(allele1);
         alleles.Push(allele2);
         }
   };

class SetOfUnknowns : public StringMap
   {
   public:
	   TUnknown * GetUnknown(int index);
      MissingAllele * GetMissing(const ::String & name);
      AlleleGraph * GetGraph(const ::String & name);

      ~SetOfUnknowns();
      
   private:
      static void * create_missing();
      static void * create_graph();
   };

#endif

 
