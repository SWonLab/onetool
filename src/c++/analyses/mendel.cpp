#include "analyses/fam.h"
#include "analyses/mendel.h"

namespace ONETOOL {

class founders_mean_GLs{

public:
	string snpid;
	std::map<std::string, long double> GLs;
	void set_snpid(string);
	void set_GLs(std::map<std::string, long double> GL);
};

void founders_mean_GLs::set_snpid(string snpid){
	this->snpid = snpid;
}

void founders_mean_GLs::set_GLs(std::map<std::string, long double> GL){
	this->GLs = GL;
}

// long double cart_product(Vi& in, std::unordered_map<std::string, double> Penetrance, string chrom){
// 
// 	Vd founders_GLiter;
// 	Digits het_founder;
// 	Digits hom_founder;
// 	vector<offd> offsprings_GLs;
// 	vector<offd> nonzeroTrans_offsprings_GL;
// 	Vd nonzeroTrans_offsprings_GLiter;
// 
// 	// 	Start founder iterators and store offspring GLs in a vector.
// 	for (Vi::const_iterator it = in.begin(); it != in.end(); ++it)
// 	{
// 		Digits d ={ (*it).GLs.begin(), (*it).GLs.end(), (*it).GLs.begin(), ((*it).mom=="0" && (*it).dad=="0") };
// 		offd od ={ (*it).GLs, (*it).xlg };
// 		if ((*it).mom=="0" && (*it).dad=="0" && (*it).xlg == "1")
// 			hom_founder = d;
// 		else if ((*it).mom=="0" && (*it).dad=="0" && (*it).xlg == "2")
// 			het_founder = d;
// 		else if ((*it).mom=="0" && (*it).dad=="0" && (*it).xlg == "0")
// 			cout << "Found parent with unknown gametic type" << endl;
// 		else
// 			offsprings_GLs.push_back(od);
// 	}
// 	founders_GLiter.push_back(hom_founder);
// 	founders_GLiter.push_back(het_founder);
// 
// 	long double Likelihood = 0.0;
// 	int fnder = 0;
// 	while (1)
// 	{
// 
// 		long double founder_prod = 1;
// 		string founder_genotype = "";
// 		//int run =0;
// 
// 
// 		for (Vd::const_iterator it = founders_GLiter.begin(); it != founders_GLiter.end(); it++)
// 		{
// 			founder_prod *= it->me->second;
// 			founder_genotype += it->me->first;
// 		}
// 		long double prod0 = founder_prod;
// 		for (auto noit = offsprings_GLs.begin(); noit != offsprings_GLs.end(); ++noit)
// 		{
// 			string chrom_offspring_gametic;
// 			if (chrom == "A")
// 				chrom_offspring_gametic = chrom;
// 			else if (chrom == "X")
// 				chrom_offspring_gametic = chrom+noit->xlg;
// 			long double sum = 0;
// 			for (auto it = (noit->GLs).begin(); it != (noit->GLs).end(); it++)
// 			{
// 				sum += it->second * Penetrance[founder_genotype+it->first+chrom_offspring_gametic];
// 			}
// 			prod0 *= sum;
// 		}
// 		Likelihood += prod0;
// 
// 		nonzeroTrans_offsprings_GL.clear();
// 		nonzeroTrans_offsprings_GLiter.clear();
// 		for (Vd::iterator it = founders_GLiter.end()-1; ;)
// 		{
// 			//			cout << "was here" << endl;
// 			// okay, I started at the left instead. sue me
// 			++(it->me);
// 			//			cout << "was here" << endl;
// 			if (it->me == it->end)
// 			{
// 				if (it == founders_GLiter.begin())
// 				{
// 					// I'm the last digit, and I'm about to roll
// 					//					cout << founder_genotype << "\t" << Likelihood << endl;
// 					return Likelihood;
// 				} else
// 				{
// 					// cascade
// 					it->me = it->begin;
// 					--it;
// 					fnder++;
// 					//					cout <<  founder_genotype << " : " << Likelihood <<  endl;
// 				}
// 			} else
// 			{
// 				// normal
// 				fnder++;
// 				//				cout <<  founder_genotype << " : " << Likelihood <<  endl;
// 				break;
// 			}
// 
// 		}
// 		//		cout << "was here too" << endl;
// 		//		cout << founder_genotype << "\t" << Likelihood << endl;
// 	}
// 
// 
// 
// }

/* !!! WISARD-specific analysis !!! */
#if (TOOLSET_TYPE & 0x100) == 0x100

// void new_compute_likelihood(cIO *Cp_IO, vector< vector <vector <LINE>>> snps,
// 	vector< vector< LINE>> founders,
// 	unordered_map<std::string, double> Penetrance,
// 	string filename, double alpha, string unfFLAG, string phredFLAG)
// {
// 	static int func_call = 0;
// 	func_call++;
// 
// 	std::vector<string> genotypes;
// 	genotypes.reserve(10);
// 	genotypes.insert(genotypes.end(), { "AA", "AC", "AG", "AT", "CC", "CG", "CT", "GG", "GT", "TT" });
// 
// 
// 
// 
// 	//	mean of founder alleles
// 
// 	//	cout << "Computing mean of founder genotype likelihoods" << endl;
// 
// 	//	population genotype vector from founder genotypes
// 
// 	//	cout << "Computing population genotype vector from founders" << endl;
// 
// 	std::vector<string> alleles;
// 	alleles.reserve(4);
// 	alleles.insert(alleles.end(), { "A", "C", "G", "T" });
// 
// 	vector<founders_mean_GLs> popGLs;
// 	popGLs.reserve(founders.size());
// 
// 	LINE temp_popfounder;
// 	std::map<std::string, long double> pGL;
// 
// 	std::map<std::string, long double> allelePr;
// 	for (int i=0; i<4; i++)
// 	{
// 		allelePr.insert(std::make_pair(alleles[i], 0.0));
// 	}
// 
// 
// 	std::map<std::string, long double> popGT;
// 
// 	//ofstream popGTfile (filename+"_populationGL.txt", ios::out|ios::app  );
// 	for (vector< vector<LINE> >::size_type i = 0; i < founders.size(); i++)
// 	{
// 		founders_mean_GLs popgl;
// 		int num_founders = 0;
// 		long double check_norm = 0;
// 		for (auto mpos = allelePr.begin(); mpos != allelePr.end(); mpos++)
// 			mpos->second = 0;
// 		for (vector< LINE> ::size_type j = 0; j < founders[i].size(); j++)
// 		{
// 
// 			temp_popfounder = founders[i][j];
// 
// 			//			cout << temp_founder.snpid << "\t";
// 			auto cpos = temp_popfounder.GLs.begin();
// 			if (cpos->second == -1) continue;
// 			else
// 			{
// 				pGL = phred2prob(temp_popfounder, phredFLAG);
// 				allelePr["A"] += 2*pGL["AA"] + pGL["AC"] + pGL["AG"] + pGL["AT"];
// 				allelePr["C"] += 2*pGL["CC"] + pGL["AC"] + pGL["CG"] + pGL["CT"];
// 				allelePr["G"] += 2*pGL["GG"] + pGL["AG"] + pGL["CG"] + pGL["GT"];
// 				allelePr["T"] += 2*pGL["TT"] + pGL["AT"] + pGL["CT"] + pGL["GT"];
// 				num_founders++;
// 			}
// 
// 
// 		}
// 		if (num_founders !=0)
// 		{
// 			for (auto mpos = allelePr.begin(); mpos != allelePr.end(); mpos++)
// 			{
// 				mpos->second = mpos->second/(2*num_founders);
// 				check_norm += mpos->second;
// 			}
// 		} else
// 		{
// 			for (auto mpos = allelePr.begin(); mpos != allelePr.end(); mpos++)
// 			{
// 				mpos->second = 0.1;
// 				check_norm += mpos->second;
// 			}
// 		}
// 		long double popGTtot = 0;
// 		for (auto mpos = allelePr.begin(); mpos != allelePr.end(); mpos++)
// 		{
// 			mpos->second = mpos->second/check_norm;
// 			popGTtot += mpos->second;
// 		}
// 
// 		long double prodallelePr;
// 		for (auto gtpos = genotypes.begin(); gtpos != genotypes.end(); gtpos++)
// 		{
// 			if (gtpos->substr(0, 1) == gtpos->substr(1, 1))
// 				prodallelePr = allelePr[gtpos->substr(0, 1)] * allelePr[gtpos->substr(1, 1)];
// 			else
// 				prodallelePr = 2 * allelePr[gtpos->substr(0, 1)] * allelePr[gtpos->substr(1, 1)];
// 			if (prodallelePr !=0 && unfFLAG == "false")
// 			{
// 				popGT.insert(std::make_pair(*gtpos, prodallelePr));
// 			} else if (unfFLAG == "true")
// 			{
// 				popGT.insert(std::make_pair(*gtpos, 0.1));
// 			}
// 		}
// 
// 		popgl.set_snpid(temp_popfounder.snpid);
// 		popgl.set_GLs(popGT);
// 		popGLs.emplace_back(popgl);
// 		//		cout << founders[i].size() << "\t" << fmgl.snpid << endl;
// 		//popGTfile <<  popgl.snpid << "\t" << num_founders <<  "\t" ;
// 		//		for (auto mpos = popGT.begin(); mpos != popGT.end(); mpos++)
// 		//		{
// 		//			popGTfile << mpos->first << ":" << mpos->second << ",";
// 		//		}
// 		//for (auto gt = genotypes.begin(); gt != genotypes.end(); gt++)
// 		//{
// 		//	if (popGT.find(*gt) == popGT.end())
// 		//		popGTfile << "0" << "\t";
// 		//	else
// 		//		popGTfile << popGT[*gt] << "\t";
// 		//}
// 		//popGTfile <<  popGTtot << endl;
// 		popGT.clear();
// 
// 	}
// 	//popGTfile.close();
// 	//	cout << "Done!!!" << endl;
// 
// 	//ofstream ulfile (filename+"_uninformativelikelihoods.txt", ios::out|ios::app  );
// 	vector < std::map<string, long double>> UninfLikelihoods;
// 	//	auto ufm=FMGLs.begin();
// 	auto ufm=popGLs.begin();
// 	long double ulikelihood;
// 	string ulsnpid;
// 
// 	wsUint	N_vrt	= Cp_IO->sizeVariant();
// 	mFam&	Xm_fam	= Cp_IO->getFamilyData();
// 	LOOP (i, N_vrt) {
// 		int run = 0;
// 		std::map<string, long double> fammap;
// 		FOREACH (mFam_it, Xm_fam, fam) {
// 			vector<GLPROB> glpv;
// 			int famsize = 0;
// 			int homs = 0;
// 			int hets = 0;
// 			int unks = 0;
// 			FOREACH (vInt_it, fam->second.Xv_members, ind) {
// 				GLPROB glp;
// 				glp.setelem(*ind, ufm->GLs, "UNINF", phredFLAG);
// 				glpv.emplace_back(glp);
// 				famsize++;
// 				ulsnpid = ind->snpid;
// 				if (ind->sex == "1")
// 					homs++;
// 				else if (ind->sex == "2")
// 					hets++;
// 				else if (ind->sex == "0")
// 					unks++;
// 				//				cout <<  "snpid:" << ind->snpid << " famid:" << ind->familyid << " individ:" <<ind->individualid << endl;
// 			}
// 
// 			run++;
// 			//			cout <<  run << ind->snpid << ind->familyid << ind->individualid << endl;
// 			string famkeyA = "A" + itoa(famsize) + "A";
// 			string famkeyX = "X" + itoa(famsize) + "X" + itoa(homs) + "X" + itoa(hets) + "X" + itoa(unks) + "X";
// 			if (fammap.size() != 0)
// 			{
// 				if (fammap.count(famkeyA)==0)
// 				{
// 					ulikelihood = cart_product(glpv, Penetrance, "A");
// 					fammap.insert(std::make_pair(famkeyA, ulikelihood));
// 					//ulfile << left << setw(10) << ulsnpid <<"\t"<< setw(10) << famkeyA << "\t" <<right << setw(15) << setprecision(6) << ulikelihood << endl;
// 				}
// 				if (fammap.count(famkeyX)==0)
// 				{
// 					ulikelihood = cart_product(glpv, Penetrance, "X");
// 					fammap.insert(std::make_pair(famkeyX, ulikelihood));
// 					//ulfile << left << setw(10) << "\t" << ulsnpid << setw(10) << famkeyX << "\t" << right << setw(15) << setprecision(6) << ulikelihood << endl;					
// 				}
// 			} else
// 			{
// 				ulikelihood = cart_product(glpv, Penetrance, "A");
// 				fammap.insert(std::make_pair(famkeyA, ulikelihood));
// 				//ulfile << left << setw(10) << ulsnpid <<"\t" <<  setw(10) << famkeyA <<"\t"<< right << setw(15) << setprecision(6)  << ulikelihood << endl;
// 				ulikelihood = cart_product(glpv, Penetrance, "X");
// 				fammap.insert(std::make_pair(famkeyX, ulikelihood));
// 				//ulfile << left << setw(10) << ulsnpid <<"\t"<<  setw(10) << famkeyX << "\t"<<right << setw(15) << setprecision(6)  << ulikelihood << endl; 				
// 			}
// 		}
// 		ufm++;
// 		UninfLikelihoods.emplace_back(fammap);
// 	}
// 	//ulfile.close();
// 	//	cout << "Done!!!" << endl;
// 
// 
// 
// 
// 	//	computing likelihoods  (33*(4**8)+51*(2**8)+16)
// 
// 
// 	//	cout << "Computing likelihoods" << endl;
// 
// 	//	ofstream plfile ("/home/chaitanya/Desktop/nancy/pedigreelikelihoods.txt", ios::out);
// 	//	ofstream plsfile ("/home/chaitanya/Desktop/nancy/snpScores.txt", ios::out  );
// 	//	if (func_call == 1)
// 	//	{
// 	//		ofstream plfile ("/home/chaitanya/Desktop/nancy/pedigreelikelihoods.txt", ios::out);
// 	//		plfile << "SNPID" << "\t" << "FAMID" << "\t" << "LIKELIHOOD" << "\t" << "uninf-LIKELIHOOD" << "\t" << "RATIO" << "\t" << "CPU-time" << endl;
// 	//		ofstream plsfile ("/home/chaitanya/Desktop/nancy/snpScores.txt", ios::out  );
// 	//		plsfile << "SNPID" << "\t" << "SCORE" << endl;
// 	//	}
// 	//	else
// 	//	{
// 	////		plfile.close();
// 	////		plsfile.close();
// 	ofstream plfile(filename+".pedigreelikelihoods", ios::out|ios::app);
// 	ofstream plsfile(filename+".snpScores", ios::out|ios::app);
// 	//	}
// 
// 
// 
// 	vector < vector <long double>> InfLikelihoods;
// 	auto vecp = UninfLikelihoods.begin();
// 	//	auto fm=FMGLs.begin();
// 	auto fm=popGLs.begin();
// 	long double inflik_A;
// 	long double inflik_X;
// 	string plfamid;
// 	string plsnpid;
// 	//int snp_count = 0;
// 	//long int sec;
// 	//time_t seconds;
// 	//long double pedlrtA = 0;
// 	//long double pedlrtX = 0;
// 	//long double llr= 0;
// 	//long double pval = 0;
// 
// 	//	ofstream pedlrtfile (filename+"_pedigreeLRTstatistic.txt", ios::out|ios::app );
// 
// 
// 	for (auto snp = snps.begin() ; snp != snps.end(); snp++)
// 	{
// 		auto famchrompos = (*snp).begin();
// 		auto indchrompos = (*famchrompos).begin();
// 		int pedigree = 0;
// 		long double pscoreA = 0;
// 		long double pscoreX = 0;
// 		long double pedLAut = 0;
// 		long double pedLSex = 0;
// 		//		pedlrtfile << left << setw(10) << plsnpid;
// 		for (auto fam = (*snp).begin() ; fam != (*snp).end(); fam++)
// 		{
// 			vector<GLPROB> glpv;
// 			int famsize = 0;
// 			int hets = 0;
// 			int homs = 0;
// 			int unks = 0;
// 			//sec = clock();
// 			//seconds = time (NULL);
// 			for (auto ind = (*fam).begin(); ind!=(*fam).end(); ind++)
// 			{
// 				GLPROB glp;
// 				glp.setelem(*ind, fm->GLs, "INF", phredFLAG);
// 				glpv.emplace_back(glp);
// 				famsize++;
// 				plfamid = ind->familyid;
// 				plsnpid = ind->snpid;
// 				if (ind->sex == "1")
// 					homs++;
// 				else if (ind->sex == "2")
// 					hets++;
// 				else if (ind->sex == "0")
// 					unks++;
// 			}
// 			string famkeyA = "A" + itoa(famsize) + "A";
// 			string famkeyX = "X" + itoa(famsize) + "X" + itoa(homs) + "X" + itoa(hets) + "X" + itoa(unks) + "X"	;
// 			inflik_A = cart_product(glpv, Penetrance, "A");
// 			inflik_X = cart_product(glpv, Penetrance, "X");
// 			//			plfile << left << setw(10) << plsnpid << setw(10) << plfamid << setw(15) << setprecision(6) << inflik_A << setw(15) << setprecision(6) <<  inflik_X << setw(15) << setprecision(6) << (*vecp)[famkeyA] << setw(15) << setprecision(6) << (*vecp)[famkeyX] << left << setw(15) << setprecision(6) << inflik_A/(*vecp)[famkeyA] << left << setw(15) << setprecision(6) << inflik_X/(*vecp)[famkeyX]<< right << setw(15) << setprecision(10) << ((long int) clock()-sec) << endl; //time(NULL) - seconds << endl;
// 			plfile <<  indchrompos->chromosome << "\t" <<  indchrompos->position << "\t" <<  plsnpid << "\t" <<  plfamid << "\t" <<  inflik_A << "\t" <<  inflik_X << "\t" <<  (*vecp)[famkeyA] << "\t" <<  (*vecp)[famkeyX] << "\t" << inflik_A/(*vecp)[famkeyA] << "\t" <<  inflik_X/(*vecp)[famkeyX]<< endl; //<<right << setw(15) << setprecision(10) <<  //((long int) clock()-sec) << "\n"; //time(NULL) - seconds << endl;
// 			//pedlrtA = log(inflik_A); ///(*vecp)[famkeyA]);
// 			//pedlrtX = log(inflik_X); ///(*vecp)[famkeyX]);
// 			pscoreA += log(inflik_A/(*vecp)[famkeyA]);
// 			pscoreX += log(inflik_X/(*vecp)[famkeyX]);
// 			pedLAut += log(inflik_A);
// 			pedLSex += log(inflik_X);
// 			//llr = -2*(pedlrtA-pedlrtX);
// 			//			pval = bchisqr(1, llr);
// 			//			pedlrtfile << setw(15) << setprecision(6) << pedlrtA << setw(15) << setprecision(6) << pedlrtX << setw(15) << setprecision(6) << llr << setw(15) << setprecision(6) << pval;
// 			//			pedlrtfile << left << setw(10) << plsnpid << "\t" << setw(10) <<  plfamid << "\t" << setw(15) << setprecision(6) << pedlrtA << "\t" << setw(15) << setprecision(6) << pedlrtX << "\t" << setw(15) << setprecision(6) << llr << "\t" << setw(15) << setprecision(6) << pval << "\n";
// 			pedigree++;
// 		}
// 		fm++;
// 		vecp++;
// 		//double LLRstatistic = -2*(pedLAut-pedLSex);
// 		//		double alpha = 0.067;
// 		long double normL = (alpha * exp(pedLSex) + (1-alpha) * exp(pedLAut));
// 		long double posterior_Prob_Sex_linkage = (alpha * exp(pedLSex))/normL;
// 		//		long double denomratio = alpha/(alpha+(1-alpha)*exp(pedLAut - pedLSex));
// 		//		plsfile << left  << plsnpid << right << setw(15) << setprecision(6) << pscoreA << right << setw(15) << setprecision(6) << pscoreX << right << setw(15) << setprecision(6) << LLRstatistic << right << setw(15) << setprecision(6) <<  pedigree << right << setw(15) << setprecision(6) <<  bchisqr(pedigree, LLRstatistic) << endl;
// 		//		plsfile << plsnpid << "\t" << setw(15) << setprecision(6) << pscoreA << "\t" << setw(15) << setprecision(6) << pscoreX << "\t"  << setw(15) << setprecision(6) << pedLAut << "\t"  << setw(15) << setprecision(6) << pedLSex << "\t"<< setw(15) << setprecision(6) << LLRstatistic << "\t" << setw(15) << setprecision(6) <<  pedigree << "\t" << setw(15) << setprecision(6) <<  bchisqr(pedigree, LLRstatistic) << "\n";
// 		plsfile << indchrompos->chromosome << "\t" << indchrompos->position << "\t" <<  plsnpid << "\t" << pscoreA << "\t"  << pscoreX << "\t"  <<  exp(pedLAut) << "\t"  <<  exp(pedLSex) << "\t"<<  posterior_Prob_Sex_linkage   << endl;
// 		//		cout << snp_count++ << endl;
// 		//		pedlrtfile << "\n";
// 	}
// 	plfile.close();
// 	plsfile.close();
// 	//	pedlrtfile.close();
// 	//	cout << "Done!!!" << endl;
// 
// }

cMendelAnalysis::cMendelAnalysis(cIO *Cp_inpIO, cFamStrAnalysis *Cp_inpAnaFS)
	: cAnalysis(Cp_inpIO)
{
	Cp_anaFS = Cp_inpAnaFS;
}

cMendelAnalysis::~cMendelAnalysis()
{

}

inline void val2al(xVariant& X_mk, wsUint N_val, char *S_buf)
{
	if (OPT_ENABLED(indel)) switch (N_val) {
	case 0: sprintf(S_buf, "%s	%s", X_mk.indel1, X_mk.indel1); break;
	case 1: sprintf(S_buf, "%s	%s", X_mk.indel1, X_mk.indel2); break;
	case 2: sprintf(S_buf, "%s	%s", X_mk.indel2, X_mk.indel2); break;
	default:  sprintf(S_buf, "NA	NA"); break;
	} else switch (N_val) {
	case 0: sprintf(S_buf, "%c	%c", X_mk.al1, X_mk.al1); break;
	case 1: sprintf(S_buf, "%c	%c", X_mk.al1, X_mk.al2); break;
	case 2: sprintf(S_buf, "%c	%c", X_mk.al2, X_mk.al2); break;
	default:  sprintf(S_buf, "NA	NA"); break;
	}
}

inline void val2al2(xVariant& X_mk, wsUint N_val, string& S_b1, string& S_b2)
{
	if (OPT_ENABLED(indel)) switch (N_val) {
	case 0: S_b2 = S_b1 = X_mk.indel1; break;
	case 1: S_b1 = X_mk.indel1; S_b2 = X_mk.indel2; break;
	case 2: S_b2 = S_b1 = X_mk.indel2; break;
	default: S_b2 = S_b1 = "<NA>"; break;
	} else {
		char S_v1[2] = { '\0', }, S_v2[2] = { '\0', };
		S_v1[0] = X_mk.al1;
		S_v2[0] = X_mk.al2;
		switch (N_val) {
		case 0: S_b2 = S_b1 = S_v1; break;
		case 1: S_b1 = S_v1; S_b2 = S_v2; break;
		case 2: S_b2 = S_b1 = S_v2; break;
		default: S_b2 = S_b1 = "<NA>"; break;
		}
	}
}

void _mendel(char **Na_data, string &S_FID, xNucFamily &X_nf,
	xMendelStat &X_t, cTableExporter *Cp_v, vVariant &Xv_mk)
{
	wsUint		N_mk	= (wsUint)Xv_mk.size();
	wsUint		N_pi	= X_nf.Xp_pat->N_idx;
	wsUint		N_mi	= X_nf.Xp_mat->N_idx;
	char*		Na_p	= Na_data[N_pi];
	char*		Na_m	= Na_data[N_mi];
	vSampPtr&	Xa_c	= X_nf.Xp_childs;

	/* For all 'available' childs... */
	wsUint	N_fcnt = 0, N_ffail = 0;
	FOREACH (vSampPtr_it, Xa_c, i) {
		xSample&	X_s		= *(*i);
		if (X_s.N_idx == SAMP_NODATA) continue;
		char*		Na_c	= Na_data[X_s.N_idx];
		wsUint N_cnt = 0, N_fail = 0;
		for (wsUint i=0 ; i<N_mk ; i++) {
			char	N_p	= isMissing(Na_p[i]) ? 3 : Na_p[i];
			char	N_m	= isMissing(Na_m[i]) ? 3 : Na_m[i];
			char	N_c	= Na_c[i];
			string	S_m1, S_m2, S_f1, S_f2, S_c1, S_c2;
			bool	ok	= true;

			/* Skip if both missing or child-missing */
			if ((isMissing(N_p) && isMissing(N_m)) || isMissing(N_c)) continue;

			wsUint	V	= N_p*4 + N_m;
			N_cnt++;
			X_t.Xv_marker[i]++;

			switch (V) {
/* 0 0 */		case 0:				/* 0 */
					ok = N_c == 0; break;
/* 0 1, 1 0 */	case 1: case 4:		/* 0 1 */
					ok = N_c < 2; break;
/* 0 2, 2 0 */	case 2: case 8:		/* 1 */
					ok = N_c == 1; break;
/* 0 3, 3 0 */	case 3: case 12:	/* 0 1 */
					ok = N_c != 2; break;
/* 1 1 */		case 5:				/* 0 1 2 */
					break;
/* 1 2, 2 1 */	case 6: case 9:		/* 1 2 */
					ok = N_c != 0; break;
/* 2 2 */		case 10:			/* 2 */
					ok = N_c == 2; break;
/* 2 3 3 2 */	case 11: case 14:	/* 1 2 */
					ok = N_c != 0; break;
			}
			if (!ok) {
				if (Cp_v) {
					val2al2(Xv_mk[i], N_p, S_m1, S_m2);
					val2al2(Xv_mk[i], N_m, S_f1, S_f2);
					val2al2(Xv_mk[i], N_c, S_c1, S_c2);
					Cp_v->write(10, S_FID.c_str(), X_nf.Xp_pat->S_IID.c_str(),
						X_nf.Xp_mat->S_IID.c_str(), X_s.S_IID.c_str(),
						S_m1.c_str(), S_m2.c_str(),
						S_f1.c_str(), S_f2.c_str(),
						S_c1.c_str(), S_c2.c_str());
				}
				N_fail++;
				X_t.Xv_fmarker[i]++;
			}
		}
		N_fcnt += N_cnt;
		N_ffail += N_fail;
		X_t.Xv_samp[X_s.N_idx] += N_cnt;
		X_t.Xv_fsamp[X_s.N_idx] += N_fail;
	}

	X_t.Xm_fam[S_FID] += N_fcnt;
	X_t.Xm_ffam[S_FID] += N_ffail;
}

void cMendelAnalysis::run()
{
	mFam&			Xm_f	= getIO()->getFamilyData();
	mNucFam&		Xm_nf	= Cp_anaFS->getNucFamData();
	vSampPtr&		Xv_smp	= getIO()->getSample();
	cTableExporter*	Cp_v	= NULL;
	vVariant&		Xv_mkr	= getIO()->getVariant();

	if (OPT_ENABLED(verbose))
		Cp_v = new cTableExporter("mendel.all.res", "ssssssssss",
			"Entire Mendelian error report", 0,
			10, "FID", "PATID", "MATID", "IID", "PAT_A1", "PAT_A2", "MAT_A1", "MAT_A2", "A1", "A2");
	X_t.Xv_fmarker.resize(Xv_mkr.size(), 0);
	X_t.Xv_marker.resize(Xv_mkr.size(), 0);
	X_t.Xv_fsamp.resize(Xv_smp.size(), 0);
	X_t.Xv_samp.resize(Xv_smp.size(), 0);


	/* For every nuclear family, do an analysis */
	FOREACH (mNucFam_it, Xm_nf, i) {
		xNucFamily&	X_nf = i->second;

		/* Do not anything if one of parent are missing */
		if (!X_nf.Xp_pat || !X_nf.Xp_mat || X_nf.Xp_pat->N_idx == SAMP_NODATA
			|| X_nf.Xp_mat->N_idx == SAMP_NODATA)
			continue;
		/* Get mendelian error */
		_mendel(getIO()->getGenotype(), X_nf.Xp_mat->S_FID,
			X_nf, X_t, Cp_v, Xv_mkr);
	}

	/* Make report */
	cTableExporter C_fam("mendel.family.res", "szzr", "Family-wise Mendelian error report",
		0, 4, "FID", "ERROR", "AVAIL", "PROP");
	cTableExporter C_smp("mendel.sample.res", "ssizr", "Sample-wise Mendelian error report",
		0, 5, "FID", "IID", "ERROR", "AVAIL", "PROP");
	cTableExporter C_mkr("mendel.variant.res", "iir", "Variant-wise Mendelian error report",
		1, 3, "ERROR", "AVAIL", "PROP");

	/* Family-level */
	FOREACH (mFam_it, Xm_f, i) {
		map<string,__int64>::iterator X_find = X_t.Xm_fam.find(i->second.S_FID);
		/* If no entry for this FID, error and go on */
		if (X_find == X_t.Xm_fam.end() || X_find->second == 0) {
			if (!OPT_ENABLED(remna))
				C_fam.write(4, i->second.S_FID.c_str(), WIS_I64NA, WIS_I64NA, WISARD_NAN);
			continue;
		}
		/* Do report */
		__int64 res = X_find->second;
		C_fam.write(4, i->second.S_FID.c_str(), X_t.Xm_ffam[i->second.S_FID],
			res, X_t.Xm_ffam[i->second.S_FID]/(wsReal)res);
	}

	/* Sample-level */
	FOREACH (vSampPtr_it, Xv_smp, i) {
		/* Do report */
		__int64 res = X_t.Xv_samp[(*i)->N_idx];
		if (res)
			C_smp.write(5, (*i)->S_FID.c_str(), (*i)->S_IID.c_str(),
				X_t.Xv_fsamp[(*i)->N_idx], res, X_t.Xv_fsamp[(*i)->N_idx]/(wsReal)res);
//		else if (!OPT_ENABLED(remna))
//			C_smp.write(5, (*i)->S_FID.c_str(), (*i)->S_IID.c_str(),
//				WIS_I32NA, WIS_I64NA, WISARD_NAN);
	}

	/* Marker-level */
	wsUint		j = 0;
	FOREACHDO (vVariant_it, Xv_mkr, i, j++) {
		/* Do report */
		int res = X_t.Xv_marker[j];
		if (res)
			C_mkr.writeVariant(&(*i), X_t.Xv_fmarker[j], res, X_t.Xv_fmarker[j]/(wsReal)res);
		else if (!OPT_ENABLED(remna))
			C_mkr.writeVariant(&(*i), WIS_I32NA, WIS_I32NA, WISARD_NAN);
	}
}

#endif

} // End namespace ONETOOL
