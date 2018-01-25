#pragma once
#ifndef __WISARD_PROGRAM_H__
#define __WISARD_PROGRAM_H__

#if defined(_WIN32) && defined(_DEBUG)
//#	include <vld.h>
#endif

/*
bit 0  (common feature)
bit 1  (farvat-specific feature)
bit 2  (fqls-specific feature)
bit 3  (qtest-specific feature)
bit 4  (himini-specific feature)
bit 5  (pharaoh-specific feature)
bit 6  (mfarvat-specific feature)
bit 7  (farvatx-specific feature)
bit 8  (wisard/onetool-specific feature)
bit 9  (onetool-specific feature)
bit 10 (comggi-specific feature)
*/
#define TOOLSET_WISARD	0x1ff
#define TOOLSET_FARVAT	0x003
#define TOOLSET_MFQLS	0x005
#define TOOLSET_QTEST	0x009
#define TOOLSET_HIMINI	0x011
#define TOOLSET_PHARAOH	0x031
#define TOOLSET_MFARVAT	0x051
#define TOOLSET_FARVATX	0x091
#define TOOLSET_PROMUL	0x111
#define TOOLSET_HSCGGI	0x421
#define TOOLSET_ONETOOL 0x3ff

#define FUNCTION_COMMON		0x001
#define FUNCTION_FARVAT		0x002
#define FUNCTION_MFQLS		0x004
#define FUNCTION_QTEST		0x008
#define FUNCTION_HIMINI		0x010
#define FUNCTION_PHARAOH	0x020
#define FUNCTION_MFARVAT	0x040
#define FUNCTION_FARVATX	0x080
#define FUNCTION_WISARD		0x100
#define FUNCTION_ONETOOL	0x200
#define FUNCTION_HSCGGI		0x400

#define VERSION_WISARD	"1.3.1"
#define VERSION_FARVAT	"1.0.12"
#define VERSION_MFQLS	"1.0.3"
#define VERSION_QTEST	"1.0.0"
#define VERSION_HIMINI	"1.0.1"
#define VERSION_PHARAOH	"1.0.2"
#define VERSION_MFARVAT	"1.0.2"
#define VERSION_FARVATX	"1.0.0"
#define VERSION_ONETOOL "0.1.0"
#define VERSION_PROMUL	"0.9.0"
#define VERSION_HSCGGI	"0.9.0"

/* Program type */
#define TOOLSET_TYPE	TOOLSET_WISARD

#if TOOLSET_TYPE & FUNCTION_WISARD
#	define TOOLSET_NAME		"WISARD"
#	define TOOLSET_EMAIL	"me@lsy.io"
#	define TOOLSET_VER		VERSION_WISARD
#	define TOOLSET_CITE		"<NA>"
#elif TOOLSET_TYPE & FUNCTION_HIMINI
#	define TOOLSET_NAME		"Hi-Mini"
#	define TOOLSET_EMAIL	"biznok@gmail.com"
#	define TOOLSET_VER		VERSION_HIMINI
#	define TOOLSET_CITE		"<NA>"
#elif TOOLSET_TYPE & FUNCTION_FARVAT
#	define TOOLSET_NAME		"FARVAT"
#	define TOOLSET_EMAIL	"choisk0413@gmail.com"
#	define TOOLSET_VER		VERSION_FARVAT
#	define TOOLSET_CITE		"Choi S, Lee S et al. FARVAT: a FAmily-based Rare Variant Association Test\n" \
	"Bioinformatics (2014). (PMID: 25075118)"
#elif TOOLSET_TYPE & FUNCTION_FARVATX
#	define TOOLSET_NAME		"FARVATX"
#	define TOOLSET_EMAIL	"choisk0413@gmail.com"
#	define TOOLSET_VER		VERSION_FARVATX
#	define TOOLSET_CITE		"Choi S, Lee S et al., FARVATX: Family-Based Rare Variant Association\n" \
							"  Test for X-Linked Genes, Genetic Epidemiology (2016)."
#elif TOOLSET_TYPE & FUNCTION_MFARVAT
#	define TOOLSET_NAME		"mFARVAT"
#	define TOOLSET_EMAIL	"longfeistat@gmail.com"
#	define TOOLSET_VER		VERSION_MFARVAT
#	define TOOLSET_CITE		"Wang L, Lee S et al. Family-based Rare Variant Association Analysis: \n" \
							"  a Fast and Efficient Method of Multivariate Phenotype Association Analysis\n" \
							"  Genetic Epidemiology (2016)."
#elif TOOLSET_TYPE & FUNCTION_QTEST
#	define TOOLSET_NAME		"Q-test"
#	define TOOLSET_EMAIL	"jhlee1213@gmail.com"
#	define TOOLSET_VER		VERSION_QTEST
#	define TOOLSET_CITE		"Lee J et al., Gene-set association tests for next-generation sequencing data,\n" \
							"Bioinformatics (2016)."
#elif TOOLSET_TYPE & FUNCTION_MFQLS
#	define TOOLSET_NAME		"MFQLS"
#	define TOOLSET_EMAIL	"me@lsy.io"
#	define TOOLSET_VER		VERSION_MFQLS
#	define TOOLSET_CITE		"Won S, et al. Family-based association analysis: a fast and efficient method of multivariate association analysis with multiple variants\n" \
	"BMC Bioinformatics (2015). (PMID: 25887481)"
#elif TOOLSET_TYPE & FUNCTION_PHARAOH
#	define TOOLSET_NAME		"PHARAOH"
#	define TOOLSET_EMAIL	"biznok@gmail.com"
#	define TOOLSET_VER		VERSION_PHARAOH
#	define TOOLSET_CITE		"Lee S et al., Pathway-based approach using hierarchical components \n" \
							"  of collapsed rare variants, Bioinformatics (2016)."
#elif TOOLSET_TYPE & FUNCTION_ONETOOL
#	define TOOLSET_NAME		"OneTool"
#	define TOOLSET_EMAIL	"me@lsy.io"
#	define TOOLSET_VER		VERSION_ONETOOL
#	define TOOLSET_CITE		"<NA>"
#elif TOOLSET_TYPE & FUNCTION_PROMUL
#	define TOOLSET_NAME		"PHARAOH-Multi"
#	define TOOLSET_EMAIL	"biznok@gmail.com"
#	define TOOLSET_VER		VERSION_PROMUL
#	define TOOLSET_CITE		"<NA>"
#elif TOOLSET_TYPE & FUNCTION_HSCGGI
#	define TOOLSET_NAME		"HiSCom-GGI"
#	define TOOLSET_EMAIL	"biznok@gmail.com"
#	define TOOLSET_VER		VERSION_HSCGGI
#	define TOOLSET_CITE		"<NA>"
#else
#	error "Invalid toolset defined"
#endif

#endif
