#include "global/option.h"
#include "utils/log.h"
#include "output/exporter.h"
#include "global/io.h"
#include <stdarg.h>

namespace ONETOOL {

inline wsUint entryVariantFormat(cExporter *Cp, xVariant &Xs, xTableType X_t)
{
	switch (X_t) {
	case TT_TABULAR:
		Cp->fmt("%s	%s	%d	", getChrName2(Xs.chr), Xs.name, Xs.pos);
		if (IS_ASSIGNED(annogene)) {
			if (OPT_ENABLED(indel))
				Cp->fmt("%s	%s", Xs.indel2?Xs.indel2:"<NA>", Xs.anno);
			else if (Xs.al2) Cp->fmt("%c	%s", Xs.al2, Xs.anno);
			else Cp->fmt("<NA>	%s", Xs.anno);
			return 5;
		} else {
			if (OPT_ENABLED(indel)) Cp->put(Xs.indel2?Xs.indel2:"<NA>");
			else if (Xs.al2) Cp->fmt("%c", Xs.al2);
			else Cp->put("<NA>");
			return 4;
		}
		break;
	case TT_CSV:
		Cp->fmt("\"%s\",\"%s\",%d,", getChrName2(Xs.chr), Xs.name, Xs.pos);
		if (IS_ASSIGNED(annogene)) {
			if (OPT_ENABLED(indel))
				Cp->fmt("\"%s\",\"%s\"", Xs.indel2?Xs.indel2:"<NA>", Xs.anno);
			else if (Xs.al2) Cp->fmt("\"%c\",\"%s\"", Xs.al2, Xs.anno);
			else Cp->fmt("\"<NA>\",\"%s\"", Xs.anno);
			return 5;
		} else {
			if (OPT_ENABLED(indel)) Cp->fmt("\"%s\"", Xs.indel2?Xs.indel2:"<NA>");
			else if (Xs.al2) Cp->fmt("\"%c\"", Xs.al2);
			else Cp->put("\"<NA>\"");
			return 4;
		}
		break;
	case TT_XML:
		Cp->fmt("    <" HEADER_VRT_1 ">%s</" HEADER_VRT_1 ">\n"
			"    <" HEADER_VRT_2 ">%s</" HEADER_VRT_2 ">\n"
			"    <" HEADER_VRT_3 ">%d</" HEADER_VRT_3 ">\n",
			getChrName2(Xs.chr), Xs.name, Xs.pos);
		if (IS_ASSIGNED(annogene)) {
			if (OPT_ENABLED(indel))
				Cp->fmt("    <" HEADER_VRT_4 ">%s</" HEADER_VRT_4 ">\n"
					"    <" HEADER_VRT_5 ">%s</" HEADER_VRT_5 ">\n",
					Xs.indel2?Xs.indel2:"<NA>", Xs.anno);
			else if (Xs.al2)
				Cp->fmt("    <" HEADER_VRT_4 ">%c</" HEADER_VRT_4 ">\n"
					"    <" HEADER_VRT_5 ">%s</" HEADER_VRT_5 ">\n",
					Xs.al2, Xs.anno);
			else
				Cp->fmt("    <" HEADER_VRT_4 ">NA</" HEADER_VRT_4 ">\n"
					"    <" HEADER_VRT_5 ">%s</" HEADER_VRT_5 ">\n", Xs.anno);
			return 5;
		} else {
			if (OPT_ENABLED(indel))
				Cp->fmt("    <" HEADER_VRT_4 ">%s</" HEADER_VRT_4 ">\n",
					Xs.indel2?Xs.indel2:"<NA>");
			else if (Xs.al2)
				Cp->fmt("    <" HEADER_VRT_4 ">%c</" HEADER_VRT_4 ">\n", Xs.al2);
			else
				Cp->fmt("    <" HEADER_VRT_4 ">NA</" HEADER_VRT_4 ">\n");
			return 4;
		}
		break;
	case TT_JSON:
		Cp->fmt("    '" HEADER_VRT_1 "':'%s',\n"
			"    '" HEADER_VRT_2 "':'%s',\n"
			"    '" HEADER_VRT_3 "':%d,\n",
			getChrName2(Xs.chr), Xs.name, Xs.pos);
		if (IS_ASSIGNED(annogene)) {
			if (OPT_ENABLED(indel)) {
				if (Xs.indel2)
					Cp->fmt("    '" HEADER_VRT_4 "':'%s',\n"
						"    '" HEADER_VRT_5 "':'%s',\n",
						Xs.indel2, Xs.anno);
				else
					Cp->fmt("    '" HEADER_VRT_4 "':null,\n"
					"    '" HEADER_VRT_5 "':'%s',\n", Xs.anno);
			} else if (Xs.al2)
				Cp->fmt("    '" HEADER_VRT_4 "':'%c',\n"
					"    '" HEADER_VRT_5 "':'%s',\n",
				Xs.al2, Xs.anno);
			else
				Cp->fmt("    '" HEADER_VRT_4 "':null,</" HEADER_VRT_4 ">\n"
					"    '" HEADER_VRT_5 "':'%s',\n", Xs.anno);
			return 5;
		} else {
			if (OPT_ENABLED(indel))
				Cp->fmt("    '" HEADER_VRT_4 "':'%s',\n",
					Xs.indel2?Xs.indel2:"<NA>");
			else if (Xs.al2)
				Cp->fmt("    '" HEADER_VRT_4 "':'%c',\n", Xs.al2);
			else
				Cp->fmt("    '" HEADER_VRT_4 "':null,\n");
			return 4;
		}
		break;
	}
	halt("SYSERR: Should not reach to here!");
	return 0;
}

cExporter::cExporter()
{

}

cExporter::~cExporter()
{
	//pverbose("Called me?\n");
}

cExporter* cExporter::summon(wsStrCst S_ext, char B_append/*=0*/,
	xExportType X_eType/*=ET_PLAIN*/)
{
	if (X_eType == ET_AUTO)
		X_eType = OPT_ENABLED(gz) ? ET_GZIP : ET_PLAIN;

#ifndef USE_GZ
	/* Check gzip capability */
	if (X_eType == ET_GZIP)
		halt("Can't use --gz option since this version is not supports gzip");
#endif

#ifndef USE_BGZF
	/* Check BGZF capability */
	if (X_eType == ET_BGZF) {
		LOGwarn("Can't use BGZF option since this version is not supports gzip");
		X_eType = ET_PLAIN;
	}
#endif

	switch (X_eType) {
		case ET_GZIP:
			return new cGzipExporter(S_ext, B_append, 1);
			break;
#ifdef USE_GZ
		case ET_BGZF:
			return new cBgzfExporter(S_ext, B_append, 1);
			break;
#endif
		case ET_PLAIN:
			return new cPlainExporter(S_ext, B_append, 0);
			break;
		case ET_BIN:
			return new cPlainExporter(S_ext, B_append, 1);
			break;
		default:
			halt("SYSERR: Unsupported export type [%d]\n", (int)X_eType);
			break;
	}

	/* Warning fix */
	return NULL;
}

cTableExporter::cTableExporter(wsStrCst S_ext, vStr Xv_inpHeaders,
	wsStrCst S_inpFmt, const char *S_desc/*=NULL*/, char B_variant/*=0*/)
{
	S_fmt[0]	= '\0';
	N_intPos	= 0;

	/* Checklength */
	if (strlen(S_inpFmt) != Xv_inpHeaders.size())
		halt("Header size[%d] does not match with format size[%d]",
			Xv_inpHeaders.size(), strlen(S_inpFmt));
	N_elem = (wsUint)Xv_inpHeaders.size();

	/* Add header if B_variant */
	if (B_variant) N_elem += headerVariant(Xv_headers, S_fmt);
	strcat(S_fmt, S_inpFmt);

	/* Clone headers */
	FOREACH (vStr_it, Xv_inpHeaders, i) Xv_headers.push_back(*i);
	_init(S_ext, S_desc);
}

cTableExporter::cTableExporter(wsStrCst S_ext, const char** Sa_headers,
	wsUint N_sz, wsStrCst S_inpFmt, const char *S_desc/*=NULL*/,
	char B_variant/*=0*/)
{
	S_fmt[0]	= '\0';
	N_intPos	= 0;

	/* Checklength */
	if (strlen(S_inpFmt) != N_sz)
		halt("Header size[%d] does not match with format size[%d]",
			N_sz, strlen(S_inpFmt));
	N_elem = N_sz;

	/* Add header if B_variant */
	if (B_variant) N_elem += headerVariant(Xv_headers, S_fmt);
	strcat(S_fmt, S_inpFmt);

	/* Clone headers */
	LOOP (i, N_sz) Xv_headers.push_back(Sa_headers[i]);
	_init(S_ext, S_desc);
}

cTableExporter::cTableExporter(const char *S_ext, const char *S_inpFmt,
	const char *S_desc, char B_variant, wsUint N_sz, ...)
{
	va_list H_varList;
	va_start(H_varList, N_sz);

	S_fmt[0]	= '\0';
	N_intPos	= 0;

	/* Checklength */
	if (strlen(S_inpFmt) != N_sz)
		halt("Header size[%d] does not match with format size[%d]",
			N_sz, strlen(S_inpFmt));
	N_elem = N_sz;

	/* Add header if B_variant */
	if (B_variant) N_elem += headerVariant(Xv_headers, S_fmt);
	strcat(S_fmt, S_inpFmt);

	/* Clone headers */
	LOOP (i, N_sz) Xv_headers.push_back(va_arg(H_varList, char *));
	_init(S_ext, S_desc);
}

cTableExporter::cTableExporter(const char *S_ext, const char *S_inpFmt,
	const char *S_desc, char B_variant, vStr& Sa_header)
{
	S_fmt[0]	= '\0';
	N_intPos	= 0;

	/* Checklength */
	wsUint N_sz = (wsUint)Sa_header.size();
	if (strlen(S_inpFmt) != N_sz)
		halt("Header size[%d] does not match with format size[%d]",
		N_sz, strlen(S_inpFmt));
	N_elem = N_sz;

	/* Add header if B_variant */
	if (B_variant) N_elem += headerVariant(Xv_headers, S_fmt);
	strcat(S_fmt, S_inpFmt);

	/* Clone headers */
	LOOP(i, N_sz) Xv_headers.push_back(Sa_header[i]);
	_init(S_ext, S_desc);
}

cTableExporter::~cTableExporter()
{
	if (Cp_exp) {
		switch (X_tt) { /* Add last element if necessary */
		case TT_TABULAR: /* Do nothing */
		case TT_CSV:
			break;
		case TT_XML: /* Print bottommost element only */
			Cp_exp->put("</wisardResult>\n");
			break;
		case TT_JSON: /* Print array terminator only */
			Cp_exp->put("]\n");
			break;
		}
		/* Deallocate exporter */
		delete Cp_exp;
	}
}

void cTableExporter::_init(wsStrCst S_ext, wsStrCst S_desc)
{
	/* Length check */
	if (strlen(S_fmt) != Xv_headers.size())
		halt("SYSERR: TableExporter #header[%d] != #format[%d]",
			Xv_headers.size(), strlen(S_fmt));

	/* Initialize exporter */
	Cp_exp = cExporter::summon(S_ext);

	/* Set outformat */
	wsStrCst S_outFormat = IS_ASSIGNED(outformat) ? OPT_STRING(outformat) : "tabular";
	if (!stricmp(S_outFormat, "tabular")) {
		X_tt = TT_TABULAR;
	} else if (!stricmp(S_outFormat, "json")) {
		X_tt = TT_JSON;
	} else if (!stricmp(S_outFormat, "xml")) {
		X_tt = TT_XML;
	} else if (!stricmp(S_outFormat, "csv")) {
		X_tt = TT_CSV;
	} else halt("Invalid --outformat value [%s]", S_outFormat);

	if (S_desc)
		LOGoutput("%s is exported to [%s.%s]\n", S_desc, OPT_STRING(out),
			S_ext);

	/* Add first element if necessary */
	switch (X_tt) {
	case TT_TABULAR: /* Add header only */
		LOOPi (i, (int)Xv_headers.size())
			Cp_exp->fmt((i+1)==N_elem?"%s":"%s	", Xv_headers[i].c_str());
		Cp_exp->put("\n");
		break;
	case TT_CSV: /* Add header only */
		LOOPi (i, (int)Xv_headers.size())
			Cp_exp->fmt((i+1)==N_elem?"\"%s\"":"\"%s\",", Xv_headers[i].c_str());
		Cp_exp->put("\n");
		break;
	case TT_XML: /* Print topmost element only */
		Cp_exp->put("<?xml version=\"1.0\"?>\n");
		Cp_exp->put("<wisardResult>\n");
		break;
	case TT_JSON: /* Print array starter only */
		Cp_exp->put("[\n");
		break;
	}
}

inline void _writeTabCsv(cExporter* Cp_exp, xTableType X_tt, char S_fmt,
	va_list* Hp_varList)
{
	va_list& H_varList = *Hp_varList;
	switch (S_fmt) {
	case 's': // Plain string
		if (X_tt == TT_TABULAR)
			Cp_exp->put(va_arg(H_varList, char *));
		else
			Cp_exp->fmt("\"%s\"", va_arg(H_varList, char *));
		break;
	case 'i': { // Integer
		int N_val = va_arg(H_varList, int);
		if (N_val == (int)WIS_I32NA) Cp_exp->fmt("NA");
		else Cp_exp->fmt("%d", N_val);
			  } break;
	case 'u': // Unsigned integer
		Cp_exp->fmt("%u", va_arg(H_varList, wsUint)); break;
	case 'r': { // real
		wsReal R_val = va_arg(H_varList, wsReal);
		if (NA(R_val)) Cp_exp->put("NA");
		else Cp_exp->fmt("%g", R_val);
			  } break;
	case 'z': { // __int64
		__int64 N_val = va_arg(H_varList, __int64);
		if (N_val == (__int64)WIS_I64NA) Cp_exp->put("NA");
		else Cp_exp->fmt(FMT_INT64, N_val);
			  } break;
	default: halt("TableExporter format error!");
	}
}

inline void _writeXml(cExporter* Cp_exp, wsStrCst S_col, char S_fmt,
	va_list* Hp_varList)
{
	va_list& H_varList = *Hp_varList;
	switch (S_fmt) {
	case 's': // Plain string
		Cp_exp->fmt("    <%s>%s</%s>\n", S_col,
			va_arg(H_varList, char *), S_col); break;
	case 'i': { // Integer
		int N_val = va_arg(H_varList, int);
		if (N_val == (int)WIS_I32NA) Cp_exp->fmt("    <%s></%s>\n",
			S_col, S_col);
		else Cp_exp->fmt("    <%s>%d</%s>\n", S_col, N_val, S_col);
			  } break;
	case 'u': // Unsigned integer
		Cp_exp->fmt("    <%s>%u</%s>\n", S_col,
			va_arg(H_varList, wsUint), S_col); break;
	case 'r': { // real
		wsReal R_val = va_arg(H_varList, wsReal);
		if (NA(R_val)) Cp_exp->fmt("    <%s></%s>\n", S_col, S_col);
		else Cp_exp->fmt("    <%s>%g</%s>\n", S_col, R_val, S_col);
			  } break;
	case 'z': { // int64
		__int64 N_val = va_arg(H_varList, __int64);
		if (N_val == (__int64)WIS_I64NA) Cp_exp->fmt("    <%s></%s>\n", S_col, S_col);
		else Cp_exp->fmt("    <%s>" FMT_INT64 "</%s>\n", S_col,
			N_val, S_col);
			  } break;
	default: halt("TableExporter format error!");
	}
}

inline void _writeJson(cExporter* Cp_exp, wsStrCst S_col, char S_fmt,
	va_list* Hp_varList)
{
	va_list& H_varList = *Hp_varList;
	switch (S_fmt) {
	case 's': // Plain string
		Cp_exp->fmt("    '%s':'%s',\n", S_col,
			va_arg(H_varList, char *)); break;
	case 'i': { // Integer
		int N_val = va_arg(H_varList, int);
		if (N_val == (int)WIS_I32NA) Cp_exp->fmt("    '%s':null,\n", S_col);
		else Cp_exp->fmt("    '%s':%d,\n", S_col, N_val);
			  } break;
	case 'u': // Unsigned integer
		Cp_exp->fmt("    '%s':%u,\n", S_col,
			va_arg(H_varList, wsUint)); break;
	case 'r': { // real
		wsReal R_val = va_arg(H_varList, wsReal);
		if (NA(R_val)) Cp_exp->fmt("    '%s':null,\n", S_col);
		else Cp_exp->fmt("    '%s':%g,\n", S_col, R_val);
			  } break;
	case 'z': { // Unsigned integer
		__int64 N_val = va_arg(H_varList, __int64);
		if (N_val == (__int64)WIS_I64NA) Cp_exp->fmt("	'%s':null,\n", S_col);
		else Cp_exp->fmt("    '%s':" FMT_INT64 ",\n", S_col, N_val);
			  } break;
	default: halt("TableExporter format error!");
	}
}

void cTableExporter::write(wsUint N_sz, ...)
{
	va_list	H_varList;
	va_start(H_varList, N_sz);
	if ((int)N_sz != N_elem) halt("TableExporter # element error!");

	switch (X_tt) {
	case TT_TABULAR:
	case TT_CSV:
		LOOP (i, N_sz) {
			_writeTabCsv(Cp_exp, X_tt, S_fmt[i], &H_varList);
			if ((i+1) != N_sz) Cp_exp->put("\t");
		}
		break;
	case TT_XML:
		Cp_exp->put("  <wisardRecord>\n");
		LOOP (i, N_sz)
			_writeXml(Cp_exp, Xv_headers[i].c_str(), S_fmt[i], &H_varList);
		Cp_exp->put("  </wisardRecord>");
		break;
	case TT_JSON:
		Cp_exp->put("  {\n");
		LOOP (i, N_sz)
			_writeJson(Cp_exp, Xv_headers[i].c_str(), S_fmt[i], &H_varList);
		Cp_exp->put("  }");
		break;
	default:
		halt("TableExporter exportType error!");
	}
	va_end(H_varList);
	Cp_exp->fmt("\n");
}

void cTableExporter::put(wsUint N_sz, ...)
{
	int N_isz = (int)N_sz;
	va_list	H_varList;
	va_start(H_varList, N_isz);

	/* If # remained element < N_sz */
	if ((N_elem-N_intPos) < N_isz)
		halt("Misaligned internal position");

	switch (X_tt) {
	case TT_TABULAR:
	case TT_CSV:
		LOOPi (i, N_isz) {
			_writeTabCsv(Cp_exp, X_tt, S_fmt[N_intPos++], &H_varList);
			if (N_intPos != N_elem) Cp_exp->put("\t");
		}
		break;
	case TT_XML:
		if (N_intPos == 0) Cp_exp->put("  <wisardRecord>\n");
		LOOPi (i, N_isz) {
			_writeXml(Cp_exp, Xv_headers[N_intPos].c_str(), S_fmt[N_intPos],
				&H_varList);
			N_intPos++;
		}
		if (N_intPos == N_elem) Cp_exp->put("  </wisardRecord>");
		break;
	case TT_JSON:
		if (N_intPos == 0) Cp_exp->put("  {\n");
		LOOPi (i, N_isz) {
			_writeJson(Cp_exp, Xv_headers[N_intPos].c_str(), S_fmt[N_intPos],
				&H_varList);
			N_intPos++;
		}
		if (N_intPos == N_elem) Cp_exp->put("  }");
		break;
	default:
		halt("TableExporter exportType error!");
	}
	va_end(H_varList);
	Cp_exp->fmt("\n");
}

void cTableExporter::next()
{
	if (N_elem != N_intPos)
		halt("Invalid position");
	N_intPos = 0;
}

void cTableExporter::writeVariant(xVariant* Xp_vrt, ...)
{
	va_list	H_varList;
	va_start(H_varList, Xp_vrt);

	/* Print header */
	if (X_tt == TT_XML)			Cp_exp->put("  <wisardRecord>\n");
	else if (X_tt == TT_JSON)	Cp_exp->put("  {\n");

	int i = (int)entryVariantFormat(Cp_exp, *Xp_vrt, X_tt);

	if (X_tt == TT_TABULAR) for ( ; i<N_elem ; i++) {
		switch (S_fmt[i]) {
		case 's': // Plain string
			Cp_exp->fmt("	%s", va_arg(H_varList, char *)); break;
		case 'i': { // Integer
			int N_val = va_arg(H_varList, int);
			if (N_val == (int)WIS_I32NA) Cp_exp->put("	NA");
			else Cp_exp->fmt("	%d", N_val);
			} break;
		case 'u': // Unsigned integer
			Cp_exp->fmt("	%u", va_arg(H_varList, wsUint)); break;
		case 'r': { // real
			wsReal R_val = va_arg(H_varList, wsReal);
			if (NA(R_val)) Cp_exp->put("	NA");
			else Cp_exp->fmt("	%g", R_val);
			} break;
		case 'z': { // __int64
			__int64 N_val = va_arg(H_varList, __int64);
			if (N_val == (__int64)WIS_I64NA) Cp_exp->put("	NA");
			else Cp_exp->fmt("	" FMT_INT64, N_val);
			} break;
		default: halt("TableExporter format error!");
		}
	} else if (X_tt == TT_CSV) for ( ; i<N_elem ; i++) {
		switch (S_fmt[i]) {
		case 's': // Plain string
			Cp_exp->fmt(",\"%s\"", va_arg(H_varList, char *)); break;
		case 'i': { // Integer
			int N_val = va_arg(H_varList, int);
			if (N_val == (int)WIS_I32NA) Cp_exp->put(",NA");
			else Cp_exp->fmt(",%d", N_val);
			} break;
		case 'u': // Unsigned integer
			Cp_exp->fmt(",%u", va_arg(H_varList, wsUint)); break;
		case 'r': { // real
			wsReal R_val = va_arg(H_varList, wsReal);
			if (NA(R_val)) Cp_exp->put(",\"NA\"");
			else Cp_exp->fmt(",%g", R_val);
			} break;
		case 'z': { // __int64
			__int64 N_val = va_arg(H_varList, __int64);
			if (N_val == (__int64)WIS_I64NA) Cp_exp->put(",NA");
			else Cp_exp->fmt("," FMT_INT64, N_val);
			} break;
		default: halt("TableExporter format error!");
		}
	} else if (X_tt == TT_XML) {
		for ( ; i<N_elem ; i++) {
			wsStrCst S_hdr = Xv_headers[i].c_str();
			switch (S_fmt[i]) {
			case 's': // Plain string
				Cp_exp->fmt("    <%s>%s</%s>\n", S_hdr, va_arg(H_varList, char *), S_hdr); break;
			case 'i': { // Integer
				int N_val = va_arg(H_varList, int);
				if (N_val == (int)WIS_I32NA) Cp_exp->fmt("    <%s></%s>\n", S_hdr, S_hdr);
				else Cp_exp->fmt("    <%s>%d</%s>\n", S_hdr, N_val, S_hdr);
				} break;
			case 'u': // Unsigned integer
				Cp_exp->fmt("    <%s>%u</%s>\n", S_hdr, va_arg(H_varList, wsUint), S_hdr); break;
			case 'r': { // real
				wsReal R_val = va_arg(H_varList, wsReal);
				if (NA(R_val)) Cp_exp->fmt("    <%s></%s>\n", S_hdr, S_hdr);
				else Cp_exp->fmt("    <%s>%g</%s>\n", S_hdr, R_val, S_hdr);
			} break;
			case 'z': { // __int64
				__int64 N_val = va_arg(H_varList, __int64);
				if (N_val == (__int64)WIS_I64NA) Cp_exp->fmt("    <%s></%s>\n", S_hdr, S_hdr);
				else Cp_exp->fmt("    <%s>" FMT_INT64 "</%s>\n", S_hdr, N_val, S_hdr);
				} break;
			default: halt("TableExporter format error!");
			}
		}
		Cp_exp->put("  </wisardRecord>");
	} else if (X_tt == TT_JSON) {
		for ( ; i<N_elem ; i++) {
			wsStrCst S_col = Xv_headers[i].c_str();
			switch (S_fmt[i]) {
			case 's': // Plain string
				Cp_exp->fmt("    '%s':'%s',\n", S_col,
					va_arg(H_varList, char *)); break;
			case 'i': { // Integer
				int N_val = va_arg(H_varList, int);
				if (N_val == WIS_I32NA) Cp_exp->fmt("    '%s':null,\n", S_col);
				else Cp_exp->fmt("    '%s':%d,\n", S_col, N_val);
				} break;
			case 'u': // Unsigned integer
				Cp_exp->fmt("    '%s':%u,\n", S_col,
					va_arg(H_varList, wsUint)); break;
			case 'r': { // real
				wsReal R_val = va_arg(H_varList, wsReal);
				if (NA(R_val)) Cp_exp->fmt("    '%s':null,\n", S_col);
				else Cp_exp->fmt("    '%s':%g,\n", S_col, R_val);
				} break;
			case 'z': { // __int64
				__int64 N_val = va_arg(H_varList, __int64);
				if (N_val == WIS_I64NA) Cp_exp->fmt("    '%s':null,\n", S_col);
				else Cp_exp->fmt("    '%s':" FMT_INT64 ",\n", S_col, N_val);
				} break;
			default: halt("TableExporter format error!");
			}
		}
		Cp_exp->put("  },");
	} else halt("TableExporter exportType error!");
	va_end(H_varList);
	Cp_exp->fmt("\n");
}

void cTableExporter::sampleWise(cIO *Cp_IO, wsVecCst Ra_v)
{
	vSampPtr&	Xv_ptr	= Cp_IO->getSample();
	wsUint		i		= 0;

	FOREACHDO (vSampPtr_it, Xv_ptr, I, i++) {
		/* Print out */
		write(3, (*I)->S_FID.c_str(), (*I)->S_IID.c_str(), Ra_v[i]);
	}
}

cPlainExporter::cPlainExporter(wsStrCst S_ext, char B_append/*=0*/,
	char B_bin/*=0*/)
{
	sprintf(S_fn, "%s.%s", OPT_STRING(out), S_ext);

#ifdef _DEBUG
	N_written = 0;
#endif
	if (B_bin)
		H_fp = fopen(S_fn, B_append?"ab+":"wb+");
	else
		H_fp = fopen(S_fn, B_append?"a+":"w+");

	if (H_fp == NULL)
		halt("Failed to open export file [%s]", S_fn);
}

cPlainExporter::~cPlainExporter()
{
	if (H_fp != NULL)
		fclose(H_fp);
	H_fp = NULL;
//	pverbose("[%s] closed\n", S_fn);
	S_fn[0] = '\0';
}

void* cPlainExporter::open(const char *S_ext) {
	char S_fn[512];
	if (H_fp == NULL)
		fclose(H_fp);
	sprintf(S_fn, "%s.%s", OPT_STRING(out), S_ext);
	H_fp = fopen(S_fn, "w+");
	if (H_fp == NULL)
		halt("Can't open export file [%s]", S_fn);
	setbuf(H_fp, NULL);

	return H_fp;
}

void cPlainExporter::fmt(const char *S_fmt, ...)
{
	va_list	H_varList;
	/* It will not used until this value have specific meaning */
	//	int		N_prtChrCount = 0;

	va_start(H_varList, S_fmt);
#ifdef _DEBUG
	int N_prtChrCount = vfprintf(H_fp, S_fmt, H_varList);
	if (N_prtChrCount == 0 && strlen(S_fmt))
		halt("At lest [%d] bytes should be written, but none!");
#else
	vfprintf(H_fp, S_fmt, H_varList);
#endif
	va_end(H_varList);
}

void cPlainExporter::put(const char *S_fmt)
{
	if (H_fp == NULL)
		halt_fmt(WISARD_SYST_NULL_IO);
	fputs(S_fmt, H_fp);
}

void cPlainExporter::put(const char *S_fmt, wsUint N_times)
{
	if (H_fp == NULL)
		halt_fmt(WISARD_SYST_NULL_IO);
	for (wsUint i=0 ; i<N_times ; i++)
		fputs(S_fmt, H_fp);
}

void cPlainExporter::write(const void *S_buf, wsUint N_len)
{
#ifdef _DEBUG
	N_written += N_len;
#endif
	fwrite(S_buf, sizeof(char), N_len, H_fp);
#ifdef _DEBUG
	if (ftell(H_fp) != N_written)
		halt("ERROR [%d(actual) != %d(intend)]", ftell(H_fp), N_written);
#endif
}

void cPlainExporter::sampleWise(cIO *Cp_IO, wsVecCst Ra_v)
{
	vSampPtr&	Xv_ptr	= Cp_IO->getSample();
	wsUint		i		= 0;

	FOREACHDO (vSampPtr_it, Xv_ptr, I, i++) {
		/* Print out */
		fmt("%s	%s	%g\n", (*I)->S_FID.c_str(), (*I)->S_IID.c_str(),
			Ra_v[i]);
	}
}

cGzipExporter::cGzipExporter(const char *S_ext, char B_append/*=0*/,
	char B_bin/*=0*/)
{
	sprintf(S_fn, "%s.%s", OPT_STRING(out), S_ext);

#ifndef USE_GZ
	/* Check gzip capability */
	if (OPT_ENABLED(gz)) halt("Can't use --gz option since this version is not supports gzip");
#endif

	if (B_bin)
		H_fp = myGzOpen(S_fn, B_append?"ab+":"wb+");
	else
		H_fp = myGzOpen(S_fn, B_append?"a+":"w+");

	if (H_fp == NULL)
		halt("Failed to open export file [%s]", S_fn);
}

cGzipExporter::~cGzipExporter()
{
	if (H_fp != NULL)
		myGzClose(H_fp);
	H_fp = NULL;
	pverbose("[%s] closed\n", S_fn);
	S_fn[0] = '\0';
}

void* cGzipExporter::open(const char *S_ext) {
	char S_fn[512];
	if (H_fp == NULL)
		myGzClose(H_fp);
	sprintf(S_fn, "%s.%s", OPT_STRING(out), S_ext);
	H_fp = myGzOpen(S_fn, "w+");
	if (H_fp == NULL)
		halt("Can't open export file [%s]", S_fn);

	return H_fp;
}

void cGzipExporter::fmt(const char *S_fmt, ...)
{
	va_list	H_varList;
	/* It will not used until this value have specific meaning */
	//	int		N_prtChrCount = 0;

	va_start(H_varList, S_fmt);
	/*N_prtChrCount = */myGzPrint(H_fp, S_fmt, H_varList);
	va_end(H_varList);
}

void cGzipExporter::put(const char *S_fmt)
{
	if (H_fp == NULL)
		halt_fmt(WISARD_SYST_NULL_IO);
	myGzPuts(S_fmt, H_fp);
}

void cGzipExporter::put(const char *S_fmt, wsUint N_times)
{
	if (H_fp == NULL)
		halt_fmt(WISARD_SYST_NULL_IO);
	for (wsUint i=0 ; i<N_times ; i++)
		myGzPuts(S_fmt, H_fp);
}

void cGzipExporter::sampleWise(cIO *Cp_IO, wsVecCst Ra_v)
{
	vSampPtr&	Xv_ptr	= Cp_IO->getSample();
	wsUint		i		= 0;

	FOREACHDO (vSampPtr_it, Xv_ptr, I, i++) {
		/* Print out */
		fmt("%s	%s	%g\n", (*I)->S_FID.c_str(), (*I)->S_IID.c_str(),
			Ra_v[i]);
	}
}

#ifdef USE_GZ
cBgzfExporter::cBgzfExporter(const char *S_ext, char B_append/*=0*/,
	char B_bin/*=0*/)
{
	sprintf(S_fn, "%s.%s", OPT_STRING(out), S_ext);

#ifndef USE_BGZF
	/* Check bgzf capability */
	halt("Can't use BGZF option since this version is not supports gzip");
#endif

	if (B_bin)
		H_fp = bgzf_open(S_fn, B_append ? "ab+" : "wb+");
	else
		H_fp = bgzf_open(S_fn, B_append ? "a+" : "w+");

	if (H_fp == NULL)
		halt("Failed to open export file [%s]", S_fn);
}

cBgzfExporter::~cBgzfExporter()
{
	if (H_fp != NULL)
		bgzf_close(H_fp);
	H_fp = NULL;
	pverbose("[%s] closed\n", S_fn);
	S_fn[0] = '\0';
}

void* cBgzfExporter::open(const char *S_ext) {
	char S_fn[512];
	if (H_fp == NULL)
		bgzf_close(H_fp);
	sprintf(S_fn, "%s.%s", OPT_STRING(out), S_ext);
	H_fp = bgzf_open(S_fn, "w+");
	if (H_fp == NULL)
		halt("Can't open export file [%s]", S_fn);

	return H_fp;
}

void cBgzfExporter::fmt(const char *S_fmt, ...)
{
	va_list	H_varList;
	char* S_buf = new char[8192];
	/* It will not used until this value have specific meaning */
	//	int		N_prtChrCount = 0;

	va_start(H_varList, S_fmt);
	int N_write = vsprintf(S_buf, S_fmt, H_varList);
	/*N_prtChrCount = */bgzf_write(H_fp, S_buf, N_write);
	va_end(H_varList);
}

void cBgzfExporter::put(const char *S_fmt)
{
	if (H_fp == NULL)
		halt_fmt(WISARD_SYST_NULL_IO);
	bgzf_write(H_fp, S_fmt, strlen(S_fmt)-1);
	bgzf_write(H_fp, "\n", 1);
}

void cBgzfExporter::put(const char *S_fmt, wsUint N_times)
{
	if (H_fp == NULL)
		halt_fmt(WISARD_SYST_NULL_IO);
	for (wsUint i=0 ; i<N_times ; i++) {
		bgzf_write(H_fp, S_fmt, strlen(S_fmt)-1);
		bgzf_write(H_fp, "\n", 1);
	}
}

void cBgzfExporter::sampleWise(cIO *Cp_IO, wsVecCst Ra_v)
{
	vSampPtr&	Xv_ptr	= Cp_IO->getSample();
	wsUint		i		= 0;

	FOREACHDO(vSampPtr_it, Xv_ptr, I, i++) {
		/* Print out */
		fmt("%s	%s	%g\n", (*I)->S_FID.c_str(), (*I)->S_IID.c_str(),
			Ra_v[i]);
	}
}
#endif

} // End namespace ONETOOL
