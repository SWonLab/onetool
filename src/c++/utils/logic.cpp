#include "utils/util.h"
#include "utils/logic.h"
#include <math.h>

namespace ONETOOL {

#define WS_NTYPE 6 /* Number of defined operand types */
wsStrCst Sa_otName[WS_NTYPE] = {
	"LOGICAL",
	"NUMERIC",
	"STRING",
	"VARIABLE",
	"OPERATION",
	"NA"
};

/* NA definition */
xOperand X_NA = {
	"",
	{ NULL },
	OTP_NA,
	NULL
};

wsStrCst xOperand::getType()
{
	wsUint N_v = 0x01;
	for (wsUint i=0 ; i<WS_NTYPE ; i++,N_v<<=1)
		if (X_type & N_v)
			return Sa_otName[i];
	halt("Invalid operand type [%d]", X_type);
	return NULL; /* Virtually meaningless */
}

wsStrCst xOperand::toString()
{
	switch (X_type) {
	case OTP_LOGICAL: sprintf(S_buf, B_val?"true":"false"); break;
	case OTP_NUMERIC: sprintf(S_buf, "%g", R_val); break;
	case OTP_OPERATION: sprintf(S_buf, "<OPERATION>"); break;
	case OTP_NA: sprintf(S_buf, "<N/A>"); break;
	default: sprintf(S_buf, "< ??? >"); break;
	}

	return S_buf;
}

void wsError(xOperand &X_op, int X_opt)
{
	char	S_msg[512], *Sp_prt = S_msg;
	wsUint N_val = 0x01;
	//	size_t	N_pos = 0;
	Sp_prt += sprintf(S_msg, "Right value requires ");
	for (wsUint i=0,j=0 ; i<WS_NTYPE ; i++,N_val<<=1) {
		if (X_opt & N_val) {
			if (j)
				Sp_prt += sprintf(Sp_prt, " or %s", Sa_otName[i]);
			else
				Sp_prt += sprintf(Sp_prt, "%s", Sa_otName[i]);
			j++;
		}
	}
	sprintf(Sp_prt, ", but [%s]", X_op.getType());
	halt(S_msg);
}

xOperand* fetch(xOperand *Xp_op)
{
	/* FIXME */
	return NULL;
}

inline wsStrCst getOperatorName(xOperator t)
{
	static char S_opName[64];
	switch (t) {
	case OPT_EQU:	strcpy(S_opName, "equal operator");
	case OPT_NEQU:	strcpy(S_opName, "not-equal operator");
	case OPT_GT:	strcpy(S_opName, "greater than operator");
	case OPT_GTE:	strcpy(S_opName, "greater than or equal operator");
	case OPT_LE:	strcpy(S_opName, "less than operator");
	case OPT_LEE:	strcpy(S_opName, "less than or equal operator");
	default:        strcpy(S_opName, "unknown operator");
	}
	return S_opName;
}

inline bool _compare(bool v1, bool v2, xOperator t)
{
	switch (t) {
	case OPT_EQU: return v1 == v2;
	case OPT_NEQU: return v1 != v2;
	case OPT_GT: return v1 > v2;
	case OPT_GTE: return v1 >= v2;
	case OPT_LE: return v1 < v2;
	case OPT_LEE: return v1 <= v2;
	default: halt_fmt(WISARD_INVL_OPERATOR, getOperatorName(t), "bool-bool comparison");
	}

	halt_fmt(WISARD_INVL_OPERATOR, "this action should not be happen", "comparison");
	//halt("Incompatible operator assigned");
	return false; /* Virtually meaningless */
}

inline bool _compare(bool v1, float _v2, xOperator t)
{
	bool v2 = _v2 != 0.0f;
	switch (t) {
	case OPT_EQU: return v1 == v2;
	case OPT_NEQU: return v1 != v2;
	case OPT_GT: return v1 > v2;
	case OPT_GTE: return v1 >= v2;
	case OPT_LE: return v1 < v2;
	case OPT_LEE: return v1 <= v2;
	default: halt_fmt(WISARD_INVL_OPERATOR, getOperatorName(t), "bool-float comparison");
//		halt("Incompatible operator assigned");
	}

	halt_fmt(WISARD_INVL_OPERATOR, "this action should not be happen", "comparison");
//	halt("Incompatible operator assigned");
	return false; /* Virtually meaningless */
}

inline bool _compare(float v1, float v2, xOperator t)
{
	switch (t) {
	case OPT_EQU: return v1 == v2;
	case OPT_NEQU: return v1 != v2;
	case OPT_GT: return v1 > v2;
	case OPT_GTE: return v1 >= v2;
	case OPT_LE: return v1 < v2;
	case OPT_LEE: return v1 <= v2;
	default: halt_fmt(WISARD_INVL_OPERATOR, getOperatorName(t), "float-float comparison");
//		halt("Incompatible operator assigned");
	}

	halt_fmt(WISARD_INVL_OPERATOR, "this action should not be happen", "comparison");
	//halt("Incompatible operator assigned");
	return false; /* Virtually meaningless */
}

inline bool _compare(wsStrCst v1, wsStrCst v2, xOperator t)
{
	switch (t) {
	case OPT_EQU: return stricmp(v1, v2) == 0;
	case OPT_NEQU: return stricmp(v1, v2) != 0;
	case OPT_GT: return stricmp(v1, v2) > 0;
	case OPT_GTE: return stricmp(v1, v2) >= 0;
	case OPT_LE: return stricmp(v1, v2) < 0;
	case OPT_LEE: return stricmp(v1, v2) <= 0;
	default: halt_fmt(WISARD_INVL_OPERATOR, "this action should not be happen", "string-string comparison");
//		halt("Incompatible operator assigned");
	}

	halt_fmt(WISARD_INVL_OPERATOR, "this action should not be happen", "comparison");
//	halt("Incompatible operator assigned");
	return false; /* Virtually meaningless */
}

xOperand* understand(xOperation *X_o, bool B_doNotHalt/*=false*/)
{
	xOperand *Xp_op1	= X_o->X_op1;
	xOperand *Xp_op2	= X_o->X_op2;
	xOperand *Xp_ret	= new xOperand;
	Xp_ret->V_val = NULL;

	/* Evaluate both if required */
	if (Xp_op1->X_type == OTP_OPERATION)
		Xp_op1 = understand(Xp_op1->X_val);
	else if (Xp_op1->X_type == OTP_VARIABLE)
		Xp_op1 = Xp_op1->V_val;
	if (Xp_op2) {
		if (Xp_op2->X_type == OTP_OPERATION)
			Xp_op2 = understand(Xp_op2->X_val);
		else if (Xp_op2->X_type == OTP_VARIABLE)
			Xp_op2 = Xp_op2->V_val;
	}

	/* Behave by 'operator' */
	switch (X_o->X_opt) {
	case OPT_NEGATION:
		/* op2 must be NULL */
		if (Xp_op2 != NULL) halt("Right-value have two operators for NEGATION");
		/* op1 must be LOGICAL */
		ASSERT_OPERAND(*Xp_op1, OTP_LOGICAL);

		Xp_ret->X_type	= OTP_LOGICAL;
		Xp_ret->B_val	= !Xp_op1->B_val;
		break;
	case OPT_AND: case OPT_OR: {
		/* op1 must be LOGICAL or NUMERIC */
		ASSERT_OPERAND(*Xp_op1, OTP_LOGICAL|OTP_NUMERIC);
		/* op2 must not be NULL and be LOGICAL or NUMERIC */
		if (Xp_op2 == NULL) halt("Second operator of right-value is NULL");
		ASSERT_OPERAND(*Xp_op2, OTP_LOGICAL|OTP_NUMERIC);

		float r1	= Xp_op1->X_type==OTP_LOGICAL ? (float)Xp_op1->B_val : Xp_op1->R_val;
		float r2	= Xp_op2->X_type==OTP_LOGICAL ? (float)Xp_op2->B_val : Xp_op2->R_val;
		Xp_ret->X_type	= OTP_LOGICAL;

		/* Do calculation */
		switch (X_o->X_opt) {
		case OPT_AND: Xp_ret->B_val = r1!=0.0f && r2!=0.0f; break;
		case OPT_OR: Xp_ret->B_val = r1!=0.0f || r2!=0.0f; break;
		default:
			halt_fmt(WISARD_INVL_OPERATOR, getOperatorName(X_o->X_opt), "and/or opreation");
//			halt("Incompatible operator assigned");
		}
	} break;
	case OPT_ADD: case OPT_SUB: case OPT_MUL: case OPT_DIV: case OPT_POWER: {
		/* op1 must be LOGICAL or NUMERIC */
		ASSERT_OPERAND(*Xp_op1, OTP_LOGICAL|OTP_NUMERIC);
		/* op2 must not be NULL and be LOGICAL or NUMERIC */
		if (Xp_op2 == NULL) halt("Second operator of right-value is NULL");
		ASSERT_OPERAND(*Xp_op2, OTP_LOGICAL|OTP_NUMERIC);

		float r1	= Xp_op1->X_type==OTP_LOGICAL ? (float)Xp_op1->B_val : Xp_op1->R_val;
		float r2	= Xp_op2->X_type==OTP_LOGICAL ? (float)Xp_op2->B_val : Xp_op2->R_val;
		Xp_ret->X_type	= OTP_NUMERIC;
		/* Do calculation */
		switch (X_o->X_opt) {
		case OPT_ADD: Xp_ret->R_val = r1 + r2; break;
		case OPT_SUB: Xp_ret->R_val = r1 - r2; break;
		case OPT_MUL: Xp_ret->R_val = r1 * r2; break;
		case OPT_DIV: Xp_ret->R_val = r1 / r2; break;
		case OPT_POWER: Xp_ret->R_val = pow(r1, r2); break;
		default:
			halt_fmt(WISARD_INVL_OPERATOR, getOperatorName(X_o->X_opt), "numerical operation");
//			halt("Incompatible operator assigned");
		}
	} break;
	case OPT_EQU: case OPT_NEQU:
	case OPT_GT: case OPT_GTE:
	case OPT_LE: case OPT_LEE:
		Xp_ret->X_type = OTP_LOGICAL;
		/* Deprecated 150203 Opt muswt be communicative, or same type */
//		if (Xp_op1->X_type != Xp_op2->X_type)
//			goto _err;
		switch (Xp_op1->X_type) {
		case OTP_LOGICAL:
			switch (Xp_op2->X_type) {
			case OTP_LOGICAL:	Xp_ret->B_val = _compare(Xp_op1->B_val, Xp_op2->B_val, X_o->X_opt); break;
			case OTP_NUMERIC:	Xp_ret->B_val = _compare(Xp_op1->B_val, Xp_op2->R_val, X_o->X_opt); break;
			case OTP_NA:		Xp_ret->X_type = OTP_NA; break;
			default: goto _err;
			} break;
		case OTP_NUMERIC:
			switch (Xp_op2->X_type) {
			case OTP_LOGICAL:	Xp_ret->B_val = _compare(Xp_op2->B_val, Xp_op1->R_val, X_o->X_opt); break;
			case OTP_NUMERIC:	Xp_ret->B_val = _compare(Xp_op1->R_val, Xp_op2->R_val, X_o->X_opt); break;
			case OTP_NA:		Xp_ret->X_type = OTP_NA; break;
			default: goto _err;
			} break;
		case OTP_STRING:
			switch (Xp_op2->X_type) {
			case OTP_STRING:	Xp_ret->B_val = _compare(Xp_op1->S_val, Xp_op2->S_val, X_o->X_opt); break;
			case OTP_NA:		Xp_ret->X_type = OTP_NA; break;
			default: goto _err;
			} break;
		case OTP_NA:
			switch (Xp_op2->X_type) {
			case OTP_NA:
				switch (X_o->X_opt) {
				case OPT_GT: case OPT_LE: case OPT_NEQU:	Xp_ret->B_val = false; break;
				case OPT_GTE: case OPT_LEE: case OPT_EQU:	Xp_ret->B_val = true; break;
				default: goto _err;
				} break;
			default:
				Xp_ret->X_type = OTP_NA;
				break;
			}
			break;
		default:
			halt("Incompatible type assigned");
		}
		/* FIXME */
		break;
	default:
		halt_fmt(WISARD_INVL_OPERATOR, getOperatorName(X_o->X_opt), "general operation");
		break;
	}

	return Xp_ret;
_err:
	halt("Incomparable two operand, between [%s] and [%s]",
		Xp_op1->getType(), Xp_op2->getType());
	return NULL; /* Virtually meaningless */
}

xOperand* str2op(char *S_p)
{
//	char *q = NULL;
	xOperand *X = new xOperand;
	str2op(S_p, *X);
	return X;
}

void str2op(char *S_p, xOperand &X)
{
	char *q = NULL;
	X.V_val = NULL;
	/* Is numeric? */
	float R_val = (float)str2dbl(S_p, &q);
	if (!q || !q[0]) {
		X.R_val	= R_val;
		X.X_type	= OTP_NUMERIC;
		return;
	}
	/* Is integer or K/M/G? */
	wsUint N_val = strtol(S_p, &q, 10);
	if (!q || !q[0]) {
		X.R_val = (float)N_val;
		X.X_type	= OTP_NUMERIC;
	} else if (!q[2] && (q[1] == 'K' || q[1] == 'k' ||
		q[1] == 'M' || q[1] == 'm' ||
		q[1] == 'G' || q[1] == 'g')) {
		X.R_val = (float)N_val;
		X.X_type	= OTP_NUMERIC;
		if (q[1] == 'K' || q[1] == 'k')
			X.R_val *= REAL_CONST(1000.0);
		else if (q[1] == 'M' || q[1] == 'm')
			X.R_val *= REAL_CONST(1000000.0);
		else if (q[1] == 'G' || q[1] == 'g')
			X.R_val *= REAL_CONST(1000000000.0);
	}
	/* is logical? */
	else if (!stricmp(S_p, "true") || !stricmp(S_p, "false")) {
		X.B_val	= S_p[0] == 't';
		X.X_type	= OTP_LOGICAL;
	} else {
		/* string */
		X.S_val	= strdup(S_p);
		X.X_type	= OTP_STRING;
	}
}

typedef enum _xParseState {
//	PS_WLEFT,
	PS_WOP,
// 	PS_WRIGHT,
// 	PS_NUMERIC,
// 	PS_END,
	PS_WTOKEN,
} xParseState;

xOperation* parse(wsStrCst Sp_expr, wsUint N_level)
{
	/* Set priority */
	map<xOperator,int>	Xm_prior;

	/* Brace */
	Xm_prior.insert(make_pair(OPT_BRACE, 0));

	/* AND/OR */
	Xm_prior.insert(make_pair(OPT_AND, 1));
	Xm_prior.insert(make_pair(OPT_OR, 1));

	/* Comparison operators */
	Xm_prior.insert(make_pair(OPT_EQU, 2));
	Xm_prior.insert(make_pair(OPT_NEQU, 2));
	Xm_prior.insert(make_pair(OPT_LE, 2));
	Xm_prior.insert(make_pair(OPT_LEE, 2));
	Xm_prior.insert(make_pair(OPT_GT, 2));
	Xm_prior.insert(make_pair(OPT_GTE, 2));

	/* + - */
	Xm_prior.insert(make_pair(OPT_ADD, 3));
	Xm_prior.insert(make_pair(OPT_SUB, 3));

	/* * / */
	Xm_prior.insert(make_pair(OPT_DIV, 4));
	Xm_prior.insert(make_pair(OPT_MUL, 4));

	/* ^ */
	Xm_prior.insert(make_pair(OPT_POWER, 5));

	vOpElem		Xv_stack;
	vOperator	Xv_opt;
	/* Copy expression */
	char *S_expr	= strdup(Sp_expr);
	xParseState	X_stat = PS_WTOKEN;

	/* Generate operation */
	//xOperation *Xp_ret = new xOperation;

	/* For each character */
	for (char *s=S_expr ; *s ; 	) {
		char *e = NULL;
		/* Skip whitespace */
		while (*s == ' ' || *s == '\t') s++;

		/* According to type */
		switch (X_stat) {
		case PS_WTOKEN:
			/* Brace? */
			if (*s == '[') {
				Xv_opt.push_back(OPT_BRACE);
				s++;
			}
			/* Number? */
			else if (*s == '-' || *s == '+' || *s == '.' || (*s >= '0' && *s <= '9')) {
				xOperand *X_op = new xOperand;
				X_op->V_val		= NULL;
				X_op->X_type	= OTP_NUMERIC;
				X_op->R_val		= (float)str2dbl(s, &e);
				/* s == e then error */
				if (s == e) halt("Unexpected character [%c] found", *s);
				xOpElem X;
				X.B_isOperand = true;
				X.Xp_op = X_op;
				Xv_stack.push_back(X);
				/* Now need operator */
				X_stat = PS_WOP;
				/* Update position */
				s = e;
			}
			/* String? */
			else if ((*s >= 'a' && *s <= 'z') || (*s >= 'A' && *s <= 'Z') || *s == '_') {
				char *S_buf = NULL;
				wsCalloc(S_buf, char, 128);
				S_buf[128 - 1] = 0x7f;

				xOperand *X_op = new xOperand;
				X_op->V_val		= NULL;
				/* Get string anyway */
				char *x=S_buf;
				while ((*s >= 'a' && *s <= 'z') || (*s >= 'A' && *s <= 'Z') || *s == '_' || (*s >='0' && *s <= '9')) {
					if (*x == 0x7f) halt("Too long string literal!");
					*(x++) = *(s++);
				}
				*x = '\0';
				/* true or false? */
				if (!stricmp(x, "true") || !stricmp(x, "false")) {
					X_op->X_type = OTP_LOGICAL;
					X_op->B_val = x[0] == 't';
					DEALLOC(S_buf);
				} else {
					X_op->X_type	= OTP_STRING;
					X_op->S_val		= S_buf;
				}
				xOpElem X;
				X.B_isOperand = true;
				X.Xp_op = X_op;
				Xv_stack.push_back(X);
				/* Now need operator */
				X_stat = PS_WOP;
			} else
				halt("Unexpected character '%c(%d)'", *s, *s);
			break;
		case PS_WOP:
			{
				xOperator	X_op = OPT_NONE;
				wsUint		N_skip = 1;

				/* Valid operator */
				switch (*s) {
				/* Closing brace? */
				case ']':
					{
						xOperator	X_op;
						bool		B_ok = false;
						while (Xv_opt.size()) {
							vOperator_rit it = Xv_opt.rbegin();
							X_op = *it;
							if (X_op == OPT_BRACE) {
								Xv_opt.erase(--(it.base()));
								B_ok = true;
								break;
							}
							xOpElem X;
							X.B_isOperand = false;
							X.X_op = X_op;
							Xv_stack.push_back(X);
							/* Erase */
							Xv_opt.erase(--(it.base()));
						}
						if (B_ok == false)
							halt("Surplus closing brace ']' found");
					} break;
				case '+': X_op = OPT_ADD; break;
				case '-': X_op = OPT_SUB; break;
				case '/': X_op = OPT_DIV; break;
				case '*': X_op = OPT_MUL; break;
				case '^': X_op = OPT_POWER; break;
				case 'a': case 'A':
					/* AND or and */
					if ((*(s+1) == 'n' || *(s+1) == 'N') &&
						(*(s+2) == 'd' || *(s+2) == 'D')) {
						X_op = OPT_AND;
						N_skip += 2;
					}
					break;
				case '&': X_op = OPT_AND; break;
				case '|': X_op = OPT_OR; break;
				case 'o': case 'O':
					/* OR or or */
					if (*(s+1) == 'r' || *(s+1) == 'R') {
						X_op = OPT_OR;
						N_skip++;
					}
					break;
				case '!':
					/* Is second character '=' */
					switch (*(s+1)) {
					case '=': X_op = OPT_NEQU; N_skip++; break;
					default: halt("Negation operator(!) is not supported");
					}
					break;
				case '=':
					/* Is second character '=' or '>' or '<' */
					switch (*(s+1)) {
					case '=': X_op = OPT_EQU; N_skip++; break;
					case '<': X_op = OPT_LEE; N_skip++; break;
					case '>': X_op = OPT_GTE; N_skip++; break;
					default: halt("Assignment operator(=) is not supported");
					}
					break;
				case '<':
					/* Is second character '=' */
					X_op = *(s+1) == '=' ? OPT_LEE : OPT_LE;
					N_skip += *(s+1) == '=';
					break;
				case '>':
					/* Is second character '=' */
					X_op = *(s+1) == '=' ? OPT_GTE : OPT_GT;
					N_skip += *(s+1) == '=';
					break;
				default:
					halt("Unexpected character '%c(%d)'", *s, *s);
				}
				/* Move to next */
				s += N_skip;
				/* Insert operators to stack until its priority is met */
				if (X_op != OPT_NONE) {
					while (Xv_opt.size()) {
						/* Get operator */
						vOperator_rit	it		= Xv_opt.rbegin();
						xOperator		X_top	= *it;

						/* Get priority */
						mOprInt_it		X_pop	= Xm_prior.find(X_op);
						mOprInt_it		X_tpop	= Xm_prior.find(X_top);
						/* Its priority is ??? */
						if (X_pop == Xm_prior.end() || X_tpop == Xm_prior.end())
							halt("Failed to find priority of operator [%d]", X_op);
						/* Stop if X_pop > X_tpop */
						if (X_pop->second > X_tpop->second) break;

						/* Insert to stack and erase */
						xOpElem X;
						X.B_isOperand = false;
						X.X_op = X_top;
						Xv_stack.push_back(X);
						Xv_opt.erase(--(it.base()));
					}
					/* Insert this token to operator stack */
					Xv_opt.push_back(X_op);
				
					X_stat = PS_WTOKEN;
				}
			} break;
		}
	}

	/* Insert remained operators by stack-order */
	RFOREACH (vOperator_rit, Xv_opt, i) {
		/* No brace allowed */
		if (*i == OPT_BRACE) halt("Unclosed brace found");
		xOpElem X;
		X.B_isOperand = false;
		X.X_op = *i;
		Xv_stack.push_back(X);
	}

	/* Now build... */
	vOpElem Xv_build;
	RFOREACH (vOpElem_rit, Xv_stack, i) {
		Xv_build.push_back(*i);
		/* GOGO */
		if (Xv_build.size() < 3)
			continue;

		/* JORIP */
		while (Xv_build.size() > 1) {
			/* STOP if last two elem are not operand */
			vOpElem_rit	it		= Xv_build.rbegin();
			xOpElem&	X_v1	= *(it++);
			xOpElem&	X_v2	= *(it++);
			if (!X_v1.B_isOperand || !X_v2.B_isOperand) break;

			/* Remove two */
			it--;
			it--;
			it = vOpElem_rit(Xv_build.erase(--(it.base())));
			it = vOpElem_rit(Xv_build.erase(--(it.base())));

			/* MUST be operator */
			xOpElem &X_v3 = *it;
			if (X_v3.B_isOperand) halt("Unexpected operand!");
			Xv_build.erase(--(it.base()));

			/* Make 'operator' operand */
			xOperation *X_o	= new xOperation;
			X_o->X_op1		= X_v1.Xp_op;
			X_o->X_op2		= X_v2.Xp_op;
			X_o->X_opt		= X_v3.X_op;

			xOperand * X_oo	= new xOperand;
			X_oo->V_val		= NULL;
			X_oo->X_type	= OTP_OPERATION;
			X_oo->X_val		= X_o;

			xOpElem X;
			X.B_isOperand	= true;
			X.Xp_op			= X_oo;
			/* Insert this 'operand' */
			Xv_build.push_back(X);
		}
	}
	Xv_build[0].Xp_op->X_val->print();

	/* DO SOMETHING */
	return Xv_build[0].Xp_op->X_val;
}

bool xOperation::parse(wsStrCst Sp_expr, wsUint N_level/*=0*/)
{
	/* Set priority */
	map<xOperator,int>	Xm_prior;

	/* Brace */
	Xm_prior.insert(make_pair(OPT_BRACE, 0));

	/* AND/OR */
	Xm_prior.insert(make_pair(OPT_AND, 1));
	Xm_prior.insert(make_pair(OPT_OR, 1));

	/* Comparison operators */
	Xm_prior.insert(make_pair(OPT_EQU, 2));
	Xm_prior.insert(make_pair(OPT_NEQU, 2));
	Xm_prior.insert(make_pair(OPT_LE, 2));
	Xm_prior.insert(make_pair(OPT_LEE, 2));
	Xm_prior.insert(make_pair(OPT_GT, 2));
	Xm_prior.insert(make_pair(OPT_GTE, 2));

	/* + - */
	Xm_prior.insert(make_pair(OPT_ADD, 3));
	Xm_prior.insert(make_pair(OPT_SUB, 3));

	/* * / */
	Xm_prior.insert(make_pair(OPT_DIV, 4));
	Xm_prior.insert(make_pair(OPT_MUL, 4));

	/* ^ */
	Xm_prior.insert(make_pair(OPT_POWER, 5));

	vector<xOpElem>		Xv_stack;
	vector<xOperator>	Xv_opt;
	/* Copy expression */
	char *S_expr	= strdup(Sp_expr);
	xParseState	X_stat = PS_WTOKEN;

	/* Generate operation */
	//xOperation *Xp_ret = new xOperation;

	/* For each character */
	for (char *s=S_expr ; *s ; 	) {
		char *e = NULL;
		/* Skip whitespace */
		while (*s == ' ' || *s == '\t') s++;

		/* According to type */
		switch (X_stat) {
		case PS_WTOKEN:
			/* Brace? */
			if (*s == '[') {
				Xv_opt.push_back(OPT_BRACE);
				s++;
			}
			/* Number? */
			else if (*s == '-' || *s == '+' || *s == '.' || (*s >= '0' && *s <= '9')) {
				xOperand *X_op = new xOperand;
				X_op->V_val		= NULL;
				X_op->X_type	= OTP_NUMERIC;
				X_op->R_val		= (float)str2dbl(s, &e);
				/* s == e then error */
				if (s == e) halt("Unexpected character [%c] found", *s);
				xOpElem X;
				X.B_isOperand = true;
				X.Xp_op = X_op;
				Xv_stack.push_back(X);
				/* Now need operator */
				X_stat = PS_WOP;
				/* Update position */
				s = e;
			}
			/* String? */
			else if ((*s >= 'a' && *s <= 'z') || (*s >= 'A' && *s <= 'Z') || *s == '_') {
				char *S_buf = NULL;
				wsCalloc(S_buf, char, 128);
				S_buf[128 - 1] = 0x7f;

				xOperand *X_op = new xOperand;
				X_op->V_val		= NULL;
				/* Get string anyway */
				char *x=S_buf;
				while ((*s >= 'a' && *s <= 'z') || (*s >= 'A' && *s <= 'Z') || *s == '_') {
					if (*x == 0x7f) halt("Too long string literal!");
					*(x++) = *(s++);
				}
				*x = '\0';
				/* true or false? */
				if (!stricmp(x, "true") || !stricmp(x, "false")) {
					X_op->X_type = OTP_LOGICAL;
					X_op->B_val = x[0] == 't';
					DEALLOC(S_buf);
				} else {
					X_op->X_type	= OTP_STRING;
					X_op->S_val		= S_buf;
				}
				xOpElem X;
				X.B_isOperand = true;
				X.Xp_op = X_op;
				Xv_stack.push_back(X);
				/* Now need operator */
				X_stat = PS_WOP;
			} else
				halt("Unexpected character '%c(%d)'", *s, *s);
			break;
		case PS_WOP:
			{
				xOperator	X_op = OPT_NONE;
				wsUint		N_skip = 1;

				/* Valid operator */
				switch (*s) {
					/* Closing brace? */
				case ']':
					{
						xOperator	X_op;
						bool		B_ok = false;
						while (Xv_opt.size()) {
							vector<xOperator>::reverse_iterator it = Xv_opt.rbegin();
							X_op = *it;
							if (X_op == OPT_BRACE) {
								Xv_opt.erase(--(it.base()));
								B_ok = true;
								break;
							}
							xOpElem X;
							X.B_isOperand = false;
							X.X_op = X_op;
							Xv_stack.push_back(X);
							/* Erase */
							Xv_opt.erase(--(it.base()));
						}
						if (B_ok == false)
							halt("Surplus closing brace ']' found");
					} break;
				case '+': X_op = OPT_ADD; break;
				case '-': X_op = OPT_SUB; break;
				case '/': X_op = OPT_DIV; break;
				case '*': X_op = OPT_MUL; break;
				case '^': X_op = OPT_POWER; break;
				case 'a': case 'A':
					/* AND or and */
					if ((*(s+1) == 'n' || *(s+1) == 'N') &&
						(*(s+2) == 'd' || *(s+2) == 'D')) {
						X_op = OPT_AND;
						N_skip += 2;
					}
					break;
				case '&': X_op = OPT_AND; break;
				case '|': X_op = OPT_OR; break;
				case 'o': case 'O':
					/* OR or or */
					if (*(s+1) == 'r' || *(s+1) == 'R') {
						X_op = OPT_OR;
						N_skip++;
					}
					break;
				case '!':
					/* Is second character '=' */
					switch (*(s+1)) {
					case '=': X_op = OPT_NEQU; N_skip++; break;
					default: halt("Negation operator(!) is not supported");
					}
					break;
				case '=':
					/* Is second character '=' or '>' or '<' */
					switch (*(s+1)) {
					case '=': X_op = OPT_EQU; N_skip++; break;
					case '<': X_op = OPT_LEE; N_skip++; break;
					case '>': X_op = OPT_GTE; N_skip++; break;
					default: halt("Assignment operator(=) is not supported");
					}
					break;
				case '<':
					/* Is second character '=' */
					X_op = *(s+1) == '=' ? OPT_LEE : OPT_LE;
					N_skip += *(s+1) == '=';
					break;
				case '>':
					/* Is second character '=' */
					X_op = *(s+1) == '=' ? OPT_GTE : OPT_GT;
					N_skip += *(s+1) == '=';
					break;
				default:
					halt("Unexpected character '%c(%d)'", *s, *s);
				}
				/* Move to next */
				s += N_skip;
				/* Insert operators to stack until its priority is met */
				if (X_op != OPT_NONE) {
					while (Xv_opt.size()) {
						/* Get operator */
						vOperator_rit	it		= Xv_opt.rbegin();
						xOperator&		X_top	= *it;
						/* Get priority */
						mOprInt_it		X_pop	= Xm_prior.find(X_op);
						mOprInt_it		X_tpop	= Xm_prior.find(X_top);
						/* Its priority is ??? */
						if (X_pop == Xm_prior.end() || X_tpop == Xm_prior.end())
							halt("Failed to find priority of operator [%d]", X_op);
						/* Stop if X_pop > X_tpop */
						if (X_pop->second > X_tpop->second) break;

						/* Insert to stack and erase */
						xOpElem X;
						X.B_isOperand	= false;
						X.X_op			= X_top;
						Xv_stack.push_back(X);
						Xv_opt.erase(--(it.base()));
					}
					/* Insert this token to operator stack */
					Xv_opt.push_back(X_op);

					X_stat = PS_WTOKEN;
				}
			} break;
		}
	}

	/* Insert remained operators by stack-order */
	RFOREACH (vOperator_rit, Xv_opt, i) {
		/* No brace allowed */
		if (*i == OPT_BRACE) halt("Unclosed brace found");
		xOpElem X;
		X.B_isOperand = false;
		X.X_op = *i;
		Xv_stack.push_back(X);
	}

	/* Now build... */
	vOpElem Xv_build;
	RFOREACH (vOpElem_rit, Xv_stack, i) {
		Xv_build.push_back(*i);
		/* GOGO */
		if (Xv_build.size() < 3)
			continue;

		/* JORIP */
		while (Xv_build.size() > 1) {
			/* STOP if last two elem are not operand */
			vOpElem_rit	it		= Xv_build.rbegin();
			xOpElem&	X_v1	= *(it++);
			xOpElem&	X_v2	= *(it++);
			if (!X_v1.B_isOperand || !X_v2.B_isOperand) break;

			/* Remove two */
			it--;
			it--;
			it = vOpElem_rit(Xv_build.erase(--(it.base())));
			it = vOpElem_rit(Xv_build.erase(--(it.base())));

			/* MUST be operator */
			xOpElem&	X_v3 = *it;
			if (X_v3.B_isOperand) halt("Unexpected operand!");
			Xv_build.erase(--(it.base()));

			/* Make 'operator' operand */
			xOperation *X_o	= new xOperation;
			X_o->X_op1		= X_v1.Xp_op;
			X_o->X_op2		= X_v2.Xp_op;
			X_o->X_opt		= X_v3.X_op;

			xOperand * X_oo	= new xOperand;
			X_oo->V_val		= NULL;
			X_oo->X_type	= OTP_OPERATION;
			X_oo->X_val		= X_o;

			xOpElem X;
			X.B_isOperand	= true;
			X.Xp_op			= X_oo;
			/* Insert this 'operand' */
			Xv_build.push_back(X);
		}
	}
	//Xv_build[0].Xp_op->X_val->print();

	X_op1	= Xv_build[0].Xp_op->X_val->X_op1;
	X_op2	= Xv_build[0].Xp_op->X_val->X_op2;
	X_opt	= Xv_build[0].Xp_op->X_val->X_opt;

	/* Remove */
	//Xv_build[0].Xp_op->X_val;

	/* DO SOMETHING */
	return true;
}

// xOperation* parse(str_c Sp_expr, wsUint N_level)
// {
// 	/* Copy expression */
// 	char *S_expr	= strdup(Sp_expr);
// 	xParseState	X_stat = PS_WLEFT;
// 
// 	/* Generate operation */
// 	xOperation *Xp_ret = new xOperation;
// 
// 	/* For each character */
// 	for (char *s=S_expr ; *s ; s++) {
// 		/* If whitespace */
// 		if (*s == ' ' || *s == '\t') {
// 			switch (X_stat) {
// 			case PS_WLEFT: /* Do nothing */
// 			case PS_WOP: /* Do nothing */
// 			case PS_WRIGHT: /* Do nothing */ break;
// 			}
// 		}
// 		/* If number */
// 		else if ((*s >= '0' && *s <= '9') || *s == '-' || *s == '+' || *s == '.') {
// 			switch (X_stat) {
// 				/* Only valid for WLEFT and WRIGHT */
// 			case PS_WLEFT:
// 			case PS_WRIGHT:
// 				char		*e		= NULL;
// 				xOperand	*Xp_rcv	= X_stat == PS_WLEFT ? Xp_ret->X_op1 : Xp_ret->X_op2;
// 				Xp_rcv->X_type	= OTP_NUMERIC;
// 				Xp_rcv->R_val	= (float)str2dbl(s, &e);
// 				X_stat = X_stat == PS_WLEFT ? PS_WOP : PS_END;
// 				s = e;
// 			}
// 		}
// 		/* If braceopen */
// 		else if (*s == '[') {
// 			switch (X_stat) {
// 				/* Only valid for WLEFT and WRIGHT */
// 			case PS_WLEFT:
// 			case PS_WRIGHT: {
// 				xOperand	*Xp_rcv = X_stat == PS_WLEFT ? Xp_ret->X_op1 : Xp_ret->X_op2;
// 				Xp_rcv->X_type	= OTP_OPERATION;
// 				Xp_rcv->X_val	= parse(s+1, N_level+1);
// 				X_stat = X_stat == PS_WLEFT ? PS_WOP : PS_END;
// 							} break;
// 			default:
// 				halt("Unexpected brace open '[' found");
// 			}
// 		}
// 		/* If braceclose */
// 		else if (*s == ']') {
// 			/* Only valid for N_level+1 */
// 			if (N_level == 0)
// 				halt("Surplus brace close ']' found");
// 			switch (X_stat) {
// 				/* Only valid when END */
// 			case PS_END: goto _ret;
// 			default:
// 				halt("Unexpected brace close ']' found");
// 			}
// 		}
// 		/* Treat as variable */
// 		else if ((*s >= 'A' && *s <= 'Z') || (*s >= 'a' && *s<= 'z') ||
// 			*s == '_') {
// 				switch (X_stat) {
// 					/* Only valid for WLEFT and WRIGHT */
// 				case PS_WLEFT:
// 				case PS_WRIGHT: {
// 					char	*S_buf = NULL;
// 					MULTI_CALLOC(S_buf, char, 256);
// 					S_buf[256 - 1] = 0x7f;
// 					char	*p		= S_buf;
// 					xOperand	*Xp_rcv = X_stat == PS_WLEFT ? Xp_ret->X_op1 : Xp_ret->X_op2;
// 
// 					/* GO until whitespace */
// 					while (*s != ' ' && *s != '\t' && ((*s >= 'A' && *s <= 'Z') ||
// 						(*s >= 'a' && *s <= 'z') || *s == '_' || *s == '-' ||
// 						(*s >= '0' && *s <= '9'))) {
// 							if (*p == 0x7f) halt("Too long literal found");
// 							*(p++) = *(s++);
// 					}
// 					s--;
// 
// 					Xp_rcv->X_type	= OTP_STRING;
// 					Xp_rcv->S_val	= S_buf;
// 					X_stat = X_stat == PS_WLEFT ? PS_WOP : PS_END;
// 								} break; 
// 				}
// 		}
// 	}
// _ret:
// 	/* Remove itself */
// 	free(S_expr);
// 	return Xp_ret;
// }

void xOperation::print(wsUint N_level/*=0*/)
{

	for (wsUint i = 0; i < N_level; i++) LOGnf("    ");
	LOG("LEFT : %s ", X_op1->getType());
	switch (X_op1->X_type) {
	case OTP_LOGICAL: LOGnf(" (%s)\n", X_op1->B_val ? "true" : "false"); break;
	case OTP_NUMERIC: LOGnf(" (%g)\n", X_op1->R_val); break;
	case OTP_OPERATION:
		LOGnf("\n");
		X_op1->X_val->print(N_level + 1);
		break;
	case OTP_STRING: LOGnf(" (%s)\n", X_op1->S_val); break;
	case OTP_VARIABLE: LOGnf(" (Variable %s)\n", X_op1->S_val); break;
	case OTP_NA: LOGnf(" <N/A>"); break;
	default: LOGnf(" (Invalid type)\n"); break;
	}
	for (wsUint i = 0; i < N_level; i++) LOGnf("    ");
	switch (X_opt) {
	case OPT_ADD: LOG("OP : +\n"); break;
	case OPT_SUB: LOG("OP : -\n"); break;
	case OPT_MUL: LOG("OP : *\n"); break;
	case OPT_DIV: LOG("OP : /\n"); break;
	case OPT_POWER: LOG("OP : ^\n"); break;
	case OPT_EQU: LOG("OP : ==\n"); break;
	case OPT_NEQU: LOG("OP : !=\n"); break;
	case OPT_GT: LOG("OP : >\n"); break;
	case OPT_GTE: LOG("OP : >=\n"); break;
	case OPT_LE: LOG("OP : <\n"); break;
	case OPT_LEE: LOG("OP : <=\n"); break;
	default: LOG("OP : Invalid\n"); break;
	}
	for (wsUint i = 0; i < N_level; i++) LOGnf("    ");
	LOG("RGHT : %s ", X_op2->getType());
	switch (X_op2->X_type) {
	case OTP_LOGICAL: LOGnf(" (%s)\n", X_op2->B_val ? "true" : "false"); break;
	case OTP_NUMERIC: LOGnf(" (%g)\n", X_op2->R_val); break;
	case OTP_OPERATION:
		LOGnf("\n");
		X_op2->X_val->print(N_level + 1);
		break;
	case OTP_STRING: LOGnf(" (%s)\n", X_op2->S_val); break;
	case OTP_VARIABLE: LOGnf(" (Variable %s)\n", X_op2->S_val); break;
	case OTP_NA: LOGnf(" <N/A>"); break;
	default: LOGnf(" (Invalid type)\n"); break;
	}
	//LOGnf("\n");
}

/* Remove operand */
void remOperand(xOperand *X_op)
{
	if (X_op == NULL || X_op == &X_NA) return;
	if (X_op->X_type == OTP_STRING || X_op->X_type == OTP_VARIABLE)
		free(X_op->S_val);
	delete X_op;
}

void xOperation::assign(wsStrCst S_varName, xOperand *Xp_val)
{
	/* if op1 == operation, do recurse */
	if (X_op1->X_type == OTP_OPERATION)
		X_op1->X_val->assign(S_varName, Xp_val);
	else if ((X_op1->X_type == OTP_STRING && !stricmp(X_op1->S_val, S_varName)) ||
		(X_op1->X_type == OTP_VARIABLE && !stricmp(X_op1->S_val, S_varName))) {
		/* Change its type to variable */
		X_op1->X_type = OTP_VARIABLE;
		/* Remove operand : i.e., this operand is either string or variable */
		remOperand(X_op1->V_val);
		/* Now V_val can use */
		X_op1->V_val = Xp_val;
	}

	/* if op2 == operation, do recurse */
	if (X_op2->X_type == OTP_OPERATION)
		X_op2->X_val->assign(S_varName, Xp_val);
	else if ((X_op2->X_type == OTP_STRING && !stricmp(X_op2->S_val, S_varName)) ||
		(X_op2->X_type == OTP_VARIABLE && !stricmp(X_op2->S_val, S_varName))) {
		/* Change its type to variable */
		X_op2->X_type = OTP_VARIABLE;
		/* Remove operand : i.e., this operand is either string or variable */
		remOperand(X_op2->V_val);
		/* Now V_val can use */
		X_op2->V_val = Xp_val;
	}
}

xOperation* xOperation::clone()
{
	xOperation* X = new xOperation;
	X->X_op1	= X_op1->clone();
	X->X_op2	= X_op2->clone();
	X->X_opt	= X_opt;

	return X;
}

xOperand* xOperand::clone()
{
	xOperand *X = new xOperand;
	X->X_type	= X_type;
	X->V_val	= NULL;
	switch (X_type) {
	case OTP_LOGICAL:	X->B_val	= B_val;			break;
	case OTP_NUMERIC:	X->R_val	= R_val;			break;
	case OTP_VARIABLE:	if (V_val) X->V_val = V_val;
	case OTP_STRING:	X->S_val	= strdup(S_val);	break;
	case OTP_OPERATION:	X->X_val	= X_val->clone();	break;
	case OTP_NA: break;
	default:
		halt("Non-clonable operand type [%d]", X_type);
	}

	return X;
}

} // End namespace ONETOOL
