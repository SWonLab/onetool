#include "global/option.h"
#include "utils/util.h"
#include "utils/genet.h"

namespace ONETOOL {

/* Corresponding complement base for given base 'A', 'C', 'G', 'T' */
char Ta_compBase[256] = {
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 0x00
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 0x10
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 0x20
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 0x30
	0, 'T', 0, 'G', 0, 0, 0, 'C', 0, 0, 0, 0, 0, 0, 'N', 0, // 0x40
	0, 0, 0, 0, 'A', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 0x54
};

/* �־��� seq�κ��� [start, end] �� �̸��� ������ ret�� ��ȯ�Ѵ�
* �޸� ret�� ������ ������ �Ҵ��� �Ǿ� �־�� ��
*/
int getPartialSeq(char *S_srcSeq, int N_startOffset, int N_endOffset, char *Sp_seqRet)
{
	if (N_startOffset < N_endOffset) {
		/* Copy sequence */
		memcpy(Sp_seqRet, S_srcSeq + N_startOffset, N_endOffset - N_startOffset + 1);
		Sp_seqRet[N_endOffset - N_startOffset + 1] = '\0';
	}
	else {
		int N_normIdx, N_revIdx;
		/* 0123
		* 3210 length=4
		* GCTA
		*/
		for (N_normIdx = N_startOffset, N_revIdx = 0; N_normIdx >= N_endOffset; N_normIdx--, N_revIdx++)
			Sp_seqRet[N_revIdx] = Ta_compBase[(wsUint)S_srcSeq[N_normIdx]];
		Sp_seqRet[N_revIdx] = '\0';
	}

	return 1;
}

int getChrAuto(char *S_buf)
{
	if (S_buf[1] && S_buf[2]) {
		LOG("Invalid chromosome definition [%s]\n", S_buf);
		exit(1);
	}

	/* For HUMAN only */
	if (S_buf[0] == 'M')
		return -1;
	else if (S_buf[0] == 'Y')
		return -1;
	else if (S_buf[0] == 'X')
		return -1;
	else {
		int tmp = atoi(S_buf);
		if (tmp == 0 || tmp > OPT_NUMBER(maxNumAutoChr))
			return -1;

		return tmp;
	}
}

int getChrAuto(int N_chr)
{
	if (N_chr == 0 || N_chr > OPT_NUMBER(maxNumAutoChr))
		return -1;

	return N_chr;
}

int getChr(char *S_buf)
{
// 150903 removed
//	static str_c S_speChrName[26] = {
//		"", "", "CH", "", "", "", "", "", "",
//		"", "", "", "MT", "", "", "", "", "",
//		"", "", "", "", "", "X", "Y", "XY"
//	};
	//	if (S_buf[1] && S_buf[2])
	//		halt("Invalid chromosome definition [%s]", S_buf);

	if ((S_buf[0] == 'c' && S_buf[1] == 'h' && S_buf[2] == 'r') ||
		(S_buf[0] == 'C' && S_buf[1] == 'H' && S_buf[2] == 'R') ||
		(S_buf[0] == 'C' && S_buf[1] == 'h' && S_buf[2] == 'r'))
		S_buf = S_buf + 3;

	/* Try to convert it integer */
	int tmp = atoi(S_buf);
	/* If unrecognized, should X, XY, Y, Mt or undefined */
	if (tmp == 0) {
		char S_fnd = '\0';
		/* XY/X? */
		if (S_buf[0] == 'X') {
			if (!S_buf[1]) S_fnd = 'X';
			else if (S_buf[1] == 'Y') S_fnd = 'Z';
		} else if (S_buf[0] == 'Y') S_fnd = 'Y';
		else if (S_buf[0] == 'M') S_fnd = 'M';
		else if (S_buf[0] == 'C') S_fnd = 'C';

		/* Return -1 if undefined */
		if (!S_fnd) return -1;
		/* Find it within nonAutoChrSeq */
		for (char *Sp=OPT_STRING(nonAutoChrSeq),i=0 ; Sp&&Sp[(int)i] ; i++)
			if (Sp[(int)i] == S_fnd)
				return (int)(NAUTO_SPECIES + i + 1);
		return -1;
	}
	/* If undefined */
	if (tmp < 0) return -1;
	/* If non-autosome */
	if (tmp > (int)NAUTO_SPECIES) {
		/* Try to get index */
		int tmp2 = tmp - (int)NAUTO_SPECIES;
		/* Exceeds range */
		if (tmp2 >= (int)strlen(OPT_STRING(nonAutoChrSeq)))
			return -1;
		/* Then fine, return 'as is' */
	}
	return tmp;


	/* For HUMAN only */
// 	if (S_buf[0] == 'M')
// 		return 26;
// 	else if (S_buf[0] == 'X' && S_buf[1] == 'Y')
// 		return 25;
// 	else if (S_buf[0] == 'Y')
// 		return 24;
// 	else if (S_buf[0] == 'X')
// 		return 23;
// 	else {
// 		int tmp = atoi(S_buf);
// 		/* If unrecognized */
// 		if (tmp == 0 || tmp > (int)NCHR_SPECIES)
// 			tmp = 0;
// 		//	halt("Invalid chromosome # definition [%s]\n", S_buf);
// 
// 		return tmp;
// 	}
}

wsStrCst getChrName2(int N_chr)
{
	static wsStrCst S_chrName[43] = {
		"Un", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11",
		"12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22",
		"23", "24", "25", "26", "27", "28", "29", "30", "31", "32",
		"33", "34", "35", "36", "37", "38", "39", "40", "41", "42" };
	static wsStrCst S_speChrName[26] = {
		"", "", "CH", "", "", "", "", "", "",
		"", "", "", "MT", "", "", "", "", "",
		"", "", "", "", "", "X", "Y", "XY"
	};
	// XYZMC
	if (N_chr < 0 || (wsUint)N_chr > NCHR_SPECIES)
		return "Un";
	else if (N_chr <= (int)NAUTO_SPECIES)
		return S_chrName[N_chr];
	int N_idx = N_chr - (int)NAUTO_SPECIES - 1;
	if (N_idx >= (int)strlen(OPT_STRING(nonAutoChrSeq)))
		return "Un";
	return S_speChrName[OPT_STRING(nonAutoChrSeq)[N_idx] - 'A'];
}

} // End namespace ONETOOL
