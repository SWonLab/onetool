#include "utils/util.h"
#include "global/io.h"
#include "utils/vis.h"
#include "global/Rconn.h"

namespace ONETOOL {

void qqVariant(wsStrCst S_inpTitle, wsStrCst S_prefix, wsMat Ra_pvals, cIO *Cp_IO)
{
	return;
#ifdef USE_R
	wsUint	N_anaSamp	= Cp_IO->sizeSample();
	wsUint	N_var		= Cp_IO->sizeVariant();
	wsUint	N_pheno		= Cp_IO->sizePheno();
	vPheno&	Xv_phe		= Cp_IO->getPhenoInfo();
	vCovar&	Xv_cov		= Cp_IO->getCovInfo();
	char S_ttl[1024];
	char S_cmd[4096], *Scp = S_cmd;
	strcpy(S_ttl, "Quantile-quantile plot");
	/* Convert path to linux style*/
#ifdef _WIN32
	char S_out[1024];
	strcpy(S_out, OPT_STRING(out));
	for (char* Sp = S_out; *Sp; Sp++)
		if (*Sp == '\\')
			*Sp = '/';
#else
	wsStrCst S_out = OPT_STRING(out);
#endif
	R().Rparse("pdf('%s.%s.qq.pdf')", S_out, S_prefix);
	for (wsUint i=0 ; i<N_pheno ; i++) {
		char *Sp = S_ttl + strlen(S_ttl);
		Sp += sprintf(Sp, " (%d samples/%d variants)\n"
			"Phenotype : %s\n", N_anaSamp, N_var, 
			Xv_phe[i].S_name.c_str());
		Sp += sprintf(Sp, "Covariate(s) :");
		if (Xv_cov.size()) FOREACH (vCovar_it, Xv_cov, P)
			Sp += sprintf(Sp, " %s", P->Sp_varName);
		else Sp += sprintf(Sp, " <NONE>");

		R().sendVector("qq", Ra_pvals[i], N_var);
		Scp += sprintf(Scp, "qqunif(qq,pch=20,main='%s');", S_ttl);
	}
	sprintf(Scp, "dev.off()");
	R().Rparse(S_cmd);
#endif
}

void hist(wsStrCst S_inpTitle, wsStrCst S_prefix, wsStrCst S_xLab, wsVec Ra_obs,
	wsUint N_obs)
{
	return;
#ifdef USE_R
	char S_ttl[1024];
	char S_cmd[4096], *Scp = S_cmd;
	sprintf(S_ttl, "Histogram of %s", S_inpTitle);
#ifdef _WIN32
	char S_out[1024];
	strcpy(S_out, OPT_STRING(out));
	for (char* Sp = S_out; *Sp; Sp++)
		if (*Sp == '\\')
			*Sp = '/';
#else
	wsStrCst S_out = OPT_STRING(out);
#endif
	R().Rparse("pdf('%s.%s.hist.pdf')", S_out, S_prefix);
	R().sendVector("obs", Ra_obs, N_obs);
	Scp += sprintf(Scp, "hist(obs,main='%s',xlab='%s');", S_ttl, S_xLab);
	sprintf(Scp, "dev.off()");
	R().Rparse(S_cmd);
#endif

}

void mhtVariant(wsStrCst S_inpTitle, wsStrCst S_prefix, wsMat Ra_pvals, int** Na_chrs, cIO* Cp_IO)
{
	return;
#ifdef USE_R
	wsUint	N_anaSamp = Cp_IO->sizeSample();
	wsUint	N_var = Cp_IO->sizeVariant();
	wsUint	N_pheno = Cp_IO->sizePheno();
	vPheno&	Xv_phe = Cp_IO->getPhenoInfo();
	vCovar&	Xv_cov = Cp_IO->getCovInfo();
	char S_ttl[1024];
	char S_cmd[4096], *Scp = S_cmd;
	strcpy(S_ttl, "Manhattan plot");
#ifdef _WIN32
	char S_out[1024];
	strcpy(S_out, OPT_STRING(out));
	for (char* Sp = S_out; *Sp; Sp++)
		if (*Sp == '\\')
			*Sp = '/';
#else
	wsStrCst S_out = OPT_STRING(out);
#endif
	R().Rparse("pdf('%s.%s.mht.pdf', width=10, height=5)", S_out, S_prefix);
//	LOOP(i, 100) LOG("%d", Na_chrs[i]);
	for (wsUint i = 0; i < N_pheno; i++) {
		char *Sp = S_ttl + strlen(S_ttl);
		Sp += sprintf(Sp, " (%d samples/%d variants)\n"
			"Phenotype : %s\n", N_anaSamp, N_var,
			Xv_phe[i].S_name.c_str());
		Sp += sprintf(Sp, "Covariate(s) :");
		if (Xv_cov.size()) FOREACH(vCovar_it, Xv_cov, P)
			Sp += sprintf(Sp, " %s", P->Sp_varName);
		else Sp += sprintf(Sp, " <NONE>");

		R().sendVector("chr", Na_chrs[i], N_var);
		R().sendVector("p1", Ra_pvals[i], N_var);
		Scp += sprintf(Scp, "mht(p1,chr,'%s');", S_ttl);
	}
	sprintf(Scp, "dev.off()");
	R().Rparse(S_cmd);
#endif
}

} // End namespace ONETOOL
