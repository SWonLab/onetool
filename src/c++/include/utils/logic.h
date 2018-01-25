#pragma once
#ifndef __WISARD_LOGIC_H__
#define __WISARD_LOGIC_H__

#include "global/common.h"
#include "utils/util.h"

namespace ONETOOL {

typedef enum _xOperandType {
	OTP_LOGICAL=0x01,	/* true / false */
	OTP_NUMERIC=0x02,	/* numerical value */
	OTP_STRING=0x04,	/* plan string */
	OTP_VARIABLE=0x08,	/* variable...? */
	OTP_OPERATION=0x10,	/* Undetermined, run-time det. */
	OTP_NA=0x20,		/* NA itself */
} xOpType;

#define ASSERT_OPERAND(d,req) if (!(((d).X_type)&(req))) wsError(d, req)

typedef enum _xOperator {
	OPT_NONE,
	OPT_NEGATION,	/* res = !o1 */
	OPT_ADD,		/* res = o1 + o2 */
	OPT_SUB,		/* res = o1 - o2 */
	OPT_MUL,		/* res = o1 * o2 */
	OPT_DIV,		/* res = o1 / o2 */
	OPT_POWER,		/* res = o1 ^ o2 */
	OPT_AND,		/* res = o1 && o2 */
	OPT_OR,			/* res = o1 || o2 */
	OPT_EQU,		/* res ?= o1 == o2 */
	OPT_NEQU,		/* res ?= o1 != o2 */
	OPT_GT,			/* res ?= o1 > o2 */
	OPT_GTE,		/* res ?= o1 >= o2 */
	OPT_LE,			/* res ?= o1 < o2 */
	OPT_LEE,		/* res ?= o1 <= o2 */
	OPT_BRACE,
} xOperator;

typedef struct _xOperation xOperation;
typedef struct _xOperand xOperand;

struct _xOperation {
	xOperand*	X_op1;
	xOperand*	X_op2;
	xOperator	X_opt;
	char*		S_text;
	void		assign(wsStrCst S_varName, xOperand *Xp_val);
	void		assignNumeric(wsStrCst S_varName, float R_val);
	void		assignString(wsStrCst S_varName, wsStrCst S_val);
	void		assignLogical(wsStrCst S_varName, bool B_val);
	xOperation*	clone();
	bool		parse(wsStrCst Sp_expr, wsUint N_level=0);
	void		print(wsUint N_level = 0);
};

typedef struct _xOpElem {
	bool		B_isOperand;
	xOperand	*Xp_op;
	xOperator	X_op;
} xOpElem;

typedef vector<xOpElem> vOpElem;
typedef vOpElem::iterator vOpElem_it;
typedef vOpElem::reverse_iterator vOpElem_rit;

struct _xOperand {
	char S_buf[64];
	union {
		bool			B_val;
		float			R_val;
		char*			S_val;
		_xOperation*	X_val;
	};
	xOpType		X_type;
	xOperand*	V_val;
	wsStrCst		getType();
	wsStrCst		toString();
	xOperand*	clone();
};

typedef map<string,xOperand>	mVar;
typedef mVar::iterator			mVar_it;

xOperation*	parse(wsStrCst Sp_expr, wsUint N_level=0);
//void		printOperation(xOperation *X_op, wsUint N_level=0);
xOperand*	understand(xOperation *X_o, bool B_doNotHalt=false);
xOperand*	str2op(char *S_p);
void		str2op(char *S_p, xOperand &X);
void		remOperand(xOperand *X_op);

typedef vector<xOperator>			vOperator;
typedef vOperator::iterator			vOperator_it;
typedef vOperator::reverse_iterator	vOperator_rit;
typedef map<xOperator,int>			mOprInt;
typedef mOprInt::iterator			mOprInt_it;

extern xOperand X_NA;

} // End namespace ONETOOL

#endif
