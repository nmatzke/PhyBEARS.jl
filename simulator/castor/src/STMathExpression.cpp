#ifndef STMATHEXPRESSION_DEF
#define STMATHEXPRESSION_DEF

#include <Rcpp.h>
#include "STMathExpression.h"

using namespace std;


#pragma mark -
#pragma mark String manipulations
#pragma mark 


// returns true if string represents a numeric value
// if string is of type <space><number><space><some stuff> then only <number> is considered
inline bool STMath_isReal(string s){
	double d;
	s.erase(s.find_last_not_of(" \f\n\t\r\v" ) + 1);
	istringstream stream(s);
	stream >> d;
	return (stream.eof() && (!stream.fail()));
}


inline double STMath_string2Double(const string &number){
	return strtod(number.c_str(), NULL);
}


template<class TYPE> 
string STMath_makeString(const TYPE &data){
	ostringstream stream;
	stream << data;
	return stream.str();
}


inline string &STMath_ltrim(std::string &haystack){
	/* old code: C++98, deprecated
	haystack.erase(haystack.begin(), std::find_if(haystack.begin(), haystack.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
	haystack.erase(haystack.begin(), std::find_if(haystack.begin(), haystack.end(), std::ptr_fun<int, int>(std::isgraph)));
	*/
	// C++11 compatible
	haystack.erase(haystack.begin(), std::find_if(haystack.begin(), haystack.end(), [](int ch) { return !std::isspace(ch); }));
	return haystack;
}

inline string &STMath_rtrim(std::string &haystack){
	/* old code: C++98, deprecated
	haystack.erase(std::find_if(haystack.rbegin(), haystack.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), haystack.end());
	haystack.erase(std::find_if(haystack.rbegin(), haystack.rend(), std::ptr_fun<int, int>(std::isgraph)).base(), haystack.end());
	*/
	// C++11 compatible
	haystack.erase(std::find_if(haystack.rbegin(), haystack.rend(), [](int ch) { return !std::isspace(ch); }).base(), haystack.end());
	return haystack;
}


inline string &STMath_trim(std::string &haystack){
    return STMath_ltrim(STMath_rtrim(haystack));
}


#pragma mark -
#pragma mark MathExpression
#pragma mark 


void MathExpression::evaluateStackEntry(long i) const{
	switch(stackFunction[i]){
	case FunctionTypeRNormal: 			stackValues[i] = random_normal(stackValues[stackArguments[i][0]], stackValues[stackArguments[i][1]]); return;
	case FunctionTypeRPoisson:			stackValues[i] = random_Poisson(stackValues[stackArguments[i][0]]); return;
	case FunctionTypeRBernoulli:		stackValues[i] = random_bernoulli(stackValues[stackArguments[i][0]]); return;
	case FunctionTypeRBinomial:			stackValues[i] = random_binomial(stackValues[stackArguments[i][0]], stackValues[stackArguments[i][1]]); return;
	case FunctionTypeRUniform:			stackValues[i] = R::runif(stackValues[stackArguments[i][0]], stackValues[stackArguments[i][1]]); return;
	case FunctionTypeRLogUniform:		stackValues[i] = random_logUniformWithinInclusive(stackValues[stackArguments[i][0]], stackValues[stackArguments[i][1]]); return;
	case FunctionTypeRCauchy:			stackValues[i] = random_Cauchy(stackValues[stackArguments[i][0]], stackValues[stackArguments[i][1]]); return;
	case FunctionTypeRChiSquared:		stackValues[i] = random_chiSquare(stackValues[stackArguments[i][0]]); return;
	case FunctionTypeRTriangular:		stackValues[i] = random_triangular(stackValues[stackArguments[i][0]], stackValues[stackArguments[i][1]], stackValues[stackArguments[i][2]]); return;
	case FunctionTypePlus: 				stackValues[i] = stackValues[stackArguments[i][0]]+stackValues[stackArguments[i][1]]; return;
	case FunctionTypeMinus:				stackValues[i] = stackValues[stackArguments[i][0]]-stackValues[stackArguments[i][1]]; return;
	case FunctionTypeMultiply:			stackValues[i] = stackValues[stackArguments[i][0]]*stackValues[stackArguments[i][1]]; return;
	case FunctionTypeDivide:			stackValues[i] = stackValues[stackArguments[i][0]]/stackValues[stackArguments[i][1]]; return;
	case FunctionTypeModulo:			stackValues[i] = fmod(stackValues[stackArguments[i][0]], stackValues[stackArguments[i][1]]); return;
	case FunctionTypeMin:				stackValues[i] = min(stackValues[stackArguments[i][0]], stackValues[stackArguments[i][1]]); return;
	case FunctionTypeMax:				stackValues[i] = max(stackValues[stackArguments[i][0]], stackValues[stackArguments[i][1]]); return;
	case FunctionTypeCos:				stackValues[i] = cos(stackValues[stackArguments[i][0]]); return;
	case FunctionTypeSin:				stackValues[i] = sin(stackValues[stackArguments[i][0]]); return;
	case FunctionTypeTan:				stackValues[i] = tan(stackValues[stackArguments[i][0]]); return;
	case FunctionTypeCot:				stackValues[i] = 1.0/tan(stackValues[stackArguments[i][0]]); return;
	case FunctionTypeAcos:				stackValues[i] = acos(stackValues[stackArguments[i][0]]); return;
	case FunctionTypeAsin:				stackValues[i] = asin(stackValues[stackArguments[i][0]]); return;
	case FunctionTypeAtan:				stackValues[i] = atan(stackValues[stackArguments[i][0]]); return;
	case FunctionTypeAtan2:				stackValues[i] = atan2(stackValues[stackArguments[i][0]], stackValues[stackArguments[i][1]]); return;
	case FunctionTypeAcot:				stackValues[i] = atan(1/stackValues[stackArguments[i][0]]); return;
	case FunctionTypeCosh:				stackValues[i] = cosh(stackValues[stackArguments[i][0]]); return;
	case FunctionTypeSinh:				stackValues[i] = sinh(stackValues[stackArguments[i][0]]); return;
	case FunctionTypeTanh:				stackValues[i] = tanh(stackValues[stackArguments[i][0]]); return;
	case FunctionTypeCoth:				stackValues[i] = 1/tanh(stackValues[stackArguments[i][0]]); return;
	case FunctionTypeExp:				stackValues[i] = exp(stackValues[stackArguments[i][0]]); return;
	case FunctionTypeLog:				stackValues[i] = log(stackValues[stackArguments[i][0]]); return;
	case FunctionTypeLog10:				stackValues[i] = log10(stackValues[stackArguments[i][0]]); return;
	case FunctionTypePow:				stackValues[i] = pow(stackValues[stackArguments[i][0]], stackValues[stackArguments[i][1]]); return;
	case FunctionTypeSqrt:				stackValues[i] = sqrt(stackValues[stackArguments[i][0]]); return;
	case FunctionTypeCeil:				stackValues[i] = ceil(stackValues[stackArguments[i][0]]); return;
	case FunctionTypeFloor:				stackValues[i] = floor(stackValues[stackArguments[i][0]]); return;
	case FunctionTypeAbs:				stackValues[i] = abs(stackValues[stackArguments[i][0]]); return;
	case FunctionTypeHeaviside:			stackValues[i] = (stackValues[stackArguments[i][0]]>0 ? 1 : 0); return;
	case FunctionTypePulse:				stackValues[i] = aux_pulse(stackValues[stackArguments[i][0]],stackValues[stackArguments[i][1]],stackValues[stackArguments[i][2]]); return;
	case FunctionTypeComb:				stackValues[i] = aux_comb(stackValues[stackArguments[i][0]], stackValues[stackArguments[i][1]], stackValues[stackArguments[i][2]], stackValues[stackArguments[i][3]]); return;
	case FunctionTypePiecewise2:		stackValues[i] = (stackValues[stackArguments[i][0]]<stackValues[stackArguments[i][1]] ? stackValues[stackArguments[i][2]] : stackValues[stackArguments[i][3]]); return;
	case FunctionTypePiecewise3:		stackValues[i] = (stackValues[stackArguments[i][0]]<stackValues[stackArguments[i][1]] ? stackValues[stackArguments[i][3]] : (stackValues[stackArguments[i][0]]<stackValues[stackArguments[i][2]] ? stackValues[stackArguments[i][4]] : stackValues[stackArguments[i][5]])); return;
	case FunctionTypeNegate:			stackValues[i] = -stackValues[stackArguments[i][0]]; return;
	case FunctionTypeEscapeNAN:			stackValues[i] = (std::isnan(stackValues[stackArguments[i][0]]) ? stackValues[stackArguments[i][1]] : stackValues[stackArguments[i][0]]); return;
	case FunctionTypeEscapeNAN2:		stackValues[i] = (std::isnan(stackValues[stackArguments[i][0]]) ? stackValues[stackArguments[i][2]] : stackValues[stackArguments[i][1]]); return;
	case FunctionTypeEscapeInf:			stackValues[i] = ((stackValues[stackArguments[i][0]]==INFTY_D) || (stackValues[stackArguments[i][0]]==-INFTY_D) ? stackValues[stackArguments[i][1]] : stackValues[stackArguments[i][0]]); return;
	case FunctionTypeEscapeInf2:		stackValues[i] = ((stackValues[stackArguments[i][0]]==INFTY_D) ? stackValues[stackArguments[i][1]] : ((stackValues[stackArguments[i][0]]==-INFTY_D) ? stackValues[stackArguments[i][2]] : stackValues[stackArguments[i][0]])); return;
	case FunctionTypeEvaluateVariable:	stackValues[i] = variableValues[stackArguments[i][0]]; return;
	case FunctionTypeNumericalConstant: return;
	case FunctionTypeUnknown: 			return;
	default: return;
	}
}



FunctionType MathExpression::str2FunctionType(const string &str){
	const string s = lowercase(str);
	if(s=="rnormal"){ return FunctionTypeRNormal; }
	else if(s=="rpoisson"){ return FunctionTypeRPoisson; }
	else if(s=="rbernoulli"){ return FunctionTypeRBernoulli; }
	else if(s=="rbinomial"){ return FunctionTypeRBinomial; }
	else if(s=="runiform"){ return FunctionTypeRUniform; }
	else if(s=="rloguniform"){ return FunctionTypeRLogUniform; }
	else if(s=="rcauchy"){ return FunctionTypeRCauchy; }
	else if(s=="rchisquared"){ return FunctionTypeRChiSquared; }
	else if(s=="rtriangular"){ return FunctionTypeRTriangular; }
	else if(s=="+"){ return FunctionTypePlus; }		
	else if(s=="-"){ return FunctionTypeMinus; }		
	else if(s=="*"){ return FunctionTypeMultiply; }	
	else if(s=="/"){ return FunctionTypeDivide; }	
	else if(s=="%"){ return FunctionTypeModulo; }	
	else if(s=="^"){ return FunctionTypePow; }		
	else if(s=="min"){ return FunctionTypeMin; }	
	else if(s=="max"){ return FunctionTypeMax; }	
	else if(s=="cos"){ return FunctionTypeCos; }		
	else if(s=="sin"){ return FunctionTypeSin; }		
	else if(s=="tan"){ return FunctionTypeTan; }		
	else if(s=="cot"){ return FunctionTypeCot; }		
	else if(s=="acos"){ return FunctionTypeAcos; }		
	else if(s=="asin"){ return FunctionTypeAsin; }		
	else if(s=="atan"){ return FunctionTypeAtan; }		
	else if(s=="atan2"){ return FunctionTypeAtan2; }		
	else if(s=="acot"){ return FunctionTypeAcot; }		
	else if(s=="cosh"){ return FunctionTypeCosh; }		
	else if(s=="sinh"){ return FunctionTypeSinh; }		
	else if(s=="tanh"){ return FunctionTypeTanh; }		
	else if(s=="coth"){ return FunctionTypeCoth; }		
	else if(s=="exp"){ return FunctionTypeExp; }		
	else if(s=="log"){ return FunctionTypeLog; }		
	else if(s=="log10"){ return FunctionTypeLog10; }		
	else if(s=="sqrt"){ return FunctionTypeSqrt; }		
	else if(s=="ceil"){ return FunctionTypeCeil; }		
	else if(s=="floor"){ return FunctionTypeFloor; }		
	else if(s=="abs"){ return FunctionTypeAbs; }		
	else if(s=="heaviside"){ return FunctionTypeHeaviside; }	
	else if(s=="pulse"){ return FunctionTypePulse; }
	else if(s=="comb"){ return FunctionTypeComb; }
	else if(s=="escapenan"){ return FunctionTypeEscapeNAN; }
	else if(s=="escapenan2"){ return FunctionTypeEscapeNAN2; }
	else if(s=="escapeinf"){ return FunctionTypeEscapeInf; }
	else if(s=="escapeinf2"){ return FunctionTypeEscapeInf2; }
	else if(s=="piecewise2"){ return FunctionTypePiecewise2; }
	else if(s=="piecewise3"){ return FunctionTypePiecewise3; }
	else{ return FunctionTypeUnknown; }
}




// returns true if expression[start:end] is whitespace (or empty)
bool MathExpression::isWhiteSpace(const string &expression, long start, long end){
	for(long i=start; i<=end; ++i){
		if(!isspace(expression[i])) return false;
	}
	return true;
}


bool MathExpression::hasWhiteSpace(const string &expression, long start, long end){
	for(long i=start; i<=end; ++i){
		if(isspace(expression[i])) return true;
	}
	return false;
}

bool MathExpression::is_one_of(char c, const char *list){
	for(long i=0; list[i]!='\0'; ++i){
		if(c==list[i]) return true;
	}
	return false;
}



// assuming pos points to a '+' or '-' in expression, follow a possible '+-' chain to the right and evaluate its compound sign
// For example, '+--+1' (assuming pos points to the first '+') will return +1
int MathExpression::getSignOfPlusMinusChain(const string &expression, long pos, long end, long &endOfChain){
	int signOfChain = 1;
	for(long i=pos; i<=end; ++i){
		if(expression[i]=='-') signOfChain = -signOfChain;
		else if(expression[i]=='+') continue;
		else{ endOfChain = i-1; return signOfChain; }
	}
	endOfChain = end;
	return signOfChain;
}



// follow plus-minus operator chains extending to the left of pos
// For example: "1 - -+5" (assuming pos points to '+') will return the position of the left-most '-'
// 				"1/  - 5" (assuming pos points to '-') will return the position of '/'
long MathExpression::getLeftMostOfPlusMinusOperatorChain(const string &expression, long start, long pos){
	long leftMost = pos;
	for(long i=pos-1; i>=start; --i){
		if(is_one_of(expression[i],"/*%^")) return i;
		else if((expression[i]=='+') || (expression[i]=='-')) leftMost = i;
		else if(!isspace(expression[i])) break;
	}
	return leftMost;
}


// search for binary operator, taking into account operator precedence (hence searching from right to left and some symbols before others)
// Elaboration: The reason to search right-to-left and in inverse operator priority is that this way the following is satisfied:
//	 			At the left-most-operator occurrence (say '+'), we can be sure that the overall expression equals <left part> + <right part> because <right part> is unsplittable at this level
// There's some outlier cases to watch out for:
//		1. 		'+' or '-' may be an exponentiation of numbers in scientific notation (e.g. "1.34e+4"). These cases will not be evaluated as operators.
// 		2a.		'+' or '-' may be "chained together" to form a "net operator" (e.g. by convention "+-56" is the same as "-56")
//		2b.		'+' or '-' may be preceded by an upper-priority operator (e.g. by convention "1/+5" is the same as "1/5")
// 		Cases (2a) and (2b) may seem excentric, but can occur in computational workflows e.g. when a symbol in an expression is blindly replaced by a negative number (e.g. "1+x" could become "1+-5")
// 		In cases (2a) & (2b) the function returns the left-most part of the chain (e.g. "+-5" returns the position of '+', and "1/+5" returns the position of '/'). This is handled by the function getLeftMostOfPlusMinusOperatorChain(..)
long MathExpression::splitBinaryOperatorInverseOrder(const string &expression, long start, long end, const std::vector<long> &pairedBrackets){
	for(long i=end; i>=start; --i){
		if(pairedBrackets[i]>=0){ i=pairedBrackets[i]; }
		else if(expression[i]=='+' || expression[i]=='-'){
			if((i<=1) || (expression[i-1]!='e' && expression[i-1]!='E')) return getLeftMostOfPlusMinusOperatorChain(expression,start,i);
			else{
				// determine if prefix in front of 'e' is a number or some other expression (e.g. variable)
				// if it is not a number, then expression[i] is indeed a binary operator, so return i
				//		otherwise it might be an exponent (scientific number notation)
				long j = expression.find_last_not_of("0123456789.",i-2);
				if(j==i-2){ return getLeftMostOfPlusMinusOperatorChain(expression,start,i); } // prefix is not a number (since 'e' is preceded by non-digit and non-dot)
				else if(j==string::npos){ if(!STMath_isReal(expression.substr(0,i-1))) return getLeftMostOfPlusMinusOperatorChain(expression,start,i); }
				else if(isValidNameCharacter(expression[j])){ return getLeftMostOfPlusMinusOperatorChain(expression,start,i); }
				else{ if(!STMath_isReal(expression.substr(j+1,i-2-j))) return getLeftMostOfPlusMinusOperatorChain(expression,start,i); }
			}
		}
	}
	for(long i=end; i>=start; --i){
		if(pairedBrackets[i]>=0){ i=pairedBrackets[i]; }
		else if(expression[i]=='*' || expression[i]=='%' || expression[i]=='/') return i;
	}
	for(long i=end; i>=start; --i){
		if(pairedBrackets[i]>=0){ i=pairedBrackets[i]; }
		else if(expression[i]=='^') return i;
	}
	return -1;
}


// returns ordered (left-to-right) list of commas separating the entries of a tuple
// comma positions are saved in commas[], and the function returns the tuple dimension (= number of commas + 1)
// expression must be a valid block
long MathExpression::splitTuple(const string &expression, long start, long end, const std::vector<long> &pairedBrackets, std::vector<long> &commas){
	commas.clear();
	for(long i=start; i<=end; ++i){
		if(pairedBrackets[i]>=0){ i=pairedBrackets[i]; }
		else if(expression[i]==','){ commas.push_back(i); }
	}
	return commas.size()+1;
}



// determines name preceding a bracket (pressumably a functio name)
// returns the left-most character of the name
// returns br if bracket has no prefix
long MathExpression::getBracketPrefix(const string &expression, long start, long br){
	if(br==start) return br;
	for(long i=br-1; i>=start; --i){
		if(!isValidNameCharacter(expression[i])) return i+1;
	}
	return start;
}








// parse an expression
// The expression must be a valid block and must be a 1-tuple
// returns an error message on failure, otherwise returns ""
string MathExpression::parseBlock(const string &expression, long start, long end, const std::vector<long> &pairedBrackets, const map<string,long> &variableName2ID, const map<string,long> &moreVariableName2ID, long stackIndex, bool allowRandoms){
	string error;
	const std::vector<long> noStackArguments;
	std::vector<long> commas;
	if(isWhiteSpace(expression,start,end)) return "Missing "+(start==0 ? "leading expression" : (end==expression.size()-1 ? "trailing expression" : "expression between '"+expression.substr(0,start)+"' and '"+expression.substr(end+1)+"'"));
		
	// find right-most binary operator at top level (if existent)
	// if such an operator (say '+') is found, then block will be equal to <left part> + <right part> because <right part> will be unsplittable at that level
	const long bo = splitBinaryOperatorInverseOrder(expression, start, end, pairedBrackets);
	
	if(bo>=0){
		// found putative binary operator at top level
		if(isWhiteSpace(expression,start,bo-1)){
			if((expression[bo]=='-') || (expression[bo]=='+')){
				long endOfChain;
				if(getSignOfPlusMinusChain(expression,bo,end,endOfChain)<0){
					// operator (or operator chain) is simply a negation of what follows
					stackFunction[stackIndex] = FunctionTypeNegate;
					stackArguments[stackIndex].push_back(stackFunction.size());
					stackFunction.push_back(FunctionTypeUnknown); // reserve spot in stack for now. FunctionType will be identified by a new call to parseBlock(..)
					stackValues.push_back(0); 
					stackArguments.push_back(noStackArguments);
					return parseBlock(expression, endOfChain+1, end, pairedBrackets, variableName2ID, moreVariableName2ID, stackFunction.size()-1, allowRandoms);
				}else if(expression[bo]=='+'){
					// operator (or operator chain) is simply a positivation of what follows
					return parseBlock(expression, endOfChain+1, end, pairedBrackets, variableName2ID, moreVariableName2ID, stackIndex, allowRandoms);
				}
			}else{
				return "Missing expression before operator '" + expression.substr(bo,1) + "'";
			}
		}
		// binary operation found at top level, so split
		stackFunction[stackIndex] = str2FunctionType(expression.substr(bo,1));
		if(stackFunction[stackIndex]==FunctionTypeUnknown) return "Unknown binary operator";
		
		stackArguments[stackIndex].push_back(stackFunction.size());
		stackFunction.push_back(FunctionTypeUnknown);
		stackValues.push_back(0); 
		stackArguments.push_back(noStackArguments);
		if((error=parseBlock(expression, start, bo-1, pairedBrackets, variableName2ID, moreVariableName2ID, stackFunction.size()-1, allowRandoms))!="") return error;
		
		stackArguments[stackIndex].push_back(stackFunction.size());
		stackFunction.push_back(FunctionTypeUnknown);
		stackValues.push_back(0); 
		stackArguments.push_back(noStackArguments);
		if((error=parseBlock(expression, bo+1, end, pairedBrackets, variableName2ID, moreVariableName2ID, stackFunction.size()-1, allowRandoms))!="") return error;

	}else{
		// no binary operation at top level
		const long rightMostBracket = expression.find_last_of(")]}", end);
		if((rightMostBracket==string::npos) || (rightMostBracket<start)){
			// this is a singleton
			string name = expression.substr(start,end-start+1);
			STMath_trim(name);
			const map<string,long>::const_iterator n2id=variableName2ID.find(name);
			const map<string,long>::const_iterator mn2id=moreVariableName2ID.find(name);			
			if((n2id!=variableName2ID.end()) || (mn2id!=moreVariableName2ID.end())){
				// expression is a valid variable
				const long ID = (n2id!=variableName2ID.end() ? n2id->second : mn2id->second);
				stackFunction[stackIndex] = FunctionTypeEvaluateVariable;
				if(ID2variable.find(ID)!=ID2variable.end()){
					// variable encountered before, so refer to it
					stackArguments[stackIndex].push_back(ID2variable[ID]);
					return "";
				}else{
					// variable encountered for the first time
					variableNames.push_back(name);
					variableIDs.push_back(ID);
					variableValues.push_back(0);
					ID2variable[ID] = variableNames.size()-1;
					stackArguments[stackIndex].push_back(variableNames.size()-1);
				}
			}else{
				// expression might be a numerical constant (or simply invalid)
				stackFunction[stackIndex] = FunctionTypeNumericalConstant;
				if(name=="pi"){ 
					stackValues[stackIndex]=3.141592653589793238462643;
				}else if(STMath_isReal(name)){
					stackValues[stackIndex]=STMath_string2Double(name);
				}else{
					return "Unknown expression '" + name + "'" + (hasWhiteSpace(name,0,name.length()-1) ? ", perhabs missing a binary operator?" : "");
				}
			}
		}else{
			// at least one bracket at this level
			if(!isWhiteSpace(expression, rightMostBracket+1, end)){
				return "Misplaced expression '" + expression.substr(rightMostBracket+1,end-rightMostBracket) + "' after bracket, perhabs missing a binary operator?";
			}
			const long pairedBracket = pairedBrackets[rightMostBracket];
			const long bp = getBracketPrefix(expression, start, pairedBracket);
			if(!isWhiteSpace(expression, start, bp-1)){
				return "Misplaced expression '" + expression.substr(start,bp-start) + "' before "+(bp==pairedBracket ? " bracket" : expression.substr(bp,pairedBracket-bp)) + ", perhabs missing a binary operator?";
			}
			if(bp==pairedBracket){
				// focal bracket is a pure bracket, so zoom in past the first bracket layer
				if((error=parseBlock(expression, pairedBracket+1, rightMostBracket-1, pairedBrackets, variableName2ID, moreVariableName2ID, stackIndex, allowRandoms))!="") return error;
				
			}else{
				// focal bracket is a function
				stackFunction[stackIndex] = str2FunctionType(expression.substr(bp,pairedBracket-bp));
				if(stackFunction[stackIndex]==FunctionTypeUnknown){
					return "Unknown function '" + expression.substr(bp,pairedBracket-bp) + "'";
				}else if((!allowRandoms) && functionTypeRandom(stackFunction[stackIndex])){
					return "Random functions such as '" + expression.substr(bp,pairedBracket-bp) + "' are not allowed";
				}
				if(functionTypeRandom(stackFunction[stackIndex])){
					hasRandoms = true;
				}
				
				// get tuple dimension (function might be multivariate)
				const long tupleDim = splitTuple(expression, pairedBracket+1, rightMostBracket-1, pairedBrackets, commas);
				if(tupleDim != functionType2dim(stackFunction[stackIndex])){
					return "Function '" + expression.substr(bp,pairedBracket-bp) + "' requires exactly " + STMath_makeString(functionType2dim(stackFunction[stackIndex])) + " arguments, but got " + STMath_makeString(tupleDim);
				}
				stackArguments[stackIndex].resize(tupleDim, -1);
				for(long c=1; c<=tupleDim; ++c){
					stackArguments[stackIndex][c-1] = stackFunction.size();
					stackFunction.push_back(FunctionTypeUnknown);
					stackValues.push_back(0); 
					stackArguments.push_back(noStackArguments);
					if((error=parseBlock(expression, (c==1 ? pairedBracket+1 : commas[c-2]+1), (c==tupleDim ? rightMostBracket-1 : commas[c-1]-1), pairedBrackets, variableName2ID, moreVariableName2ID, stackFunction.size()-1, allowRandoms))!="") return error;
				}
			}
			
		}
	}
	return "";
	
}




void MathExpression::clear(){
	stackArguments.clear();
	stackValues.clear();
	stackFunction.clear();
	variableNames.clear();
	variableValues.clear();
	variableIDs.clear();
	ID2variable.clear();
	OK = true;
	hasRandoms = false;
	originalExpression = "";
}






bool MathExpression::findBracketPairs(const string &expression, std::vector<long> &pairedBrackets){
	pairedBrackets.resize(expression.length(),-1);
	std::vector<long> openBrackets_round, openBrackets_block, openBrackets_curly;
	for(long i=0; i<expression.length(); ++i){
		if(expression[i]=='('){ openBrackets_round.push_back(i); }
		else if(expression[i]=='['){ openBrackets_block.push_back(i); }
		else if(expression[i]=='{'){ openBrackets_curly.push_back(i); }
		else if(expression[i]==')'){ 
			if(openBrackets_round.empty()) return false;
			if((!openBrackets_block.empty()) && openBrackets_block.back()>openBrackets_round.back()) return false;
			if((!openBrackets_curly.empty()) && openBrackets_curly.back()>openBrackets_round.back()) return false;
			pairedBrackets[openBrackets_round.back()] = i; 
			pairedBrackets[i] = openBrackets_round.back();
			openBrackets_round.pop_back();
		}else if(expression[i]==']'){ 
			if(openBrackets_block.empty()) return false;
			if((!openBrackets_round.empty()) && openBrackets_round.back()>openBrackets_block.back()) return false;
			if((!openBrackets_curly.empty()) && openBrackets_curly.back()>openBrackets_block.back()) return false;
			pairedBrackets[openBrackets_block.back()] = i; 
			pairedBrackets[i] = openBrackets_block.back();
			openBrackets_block.pop_back();
		}else if(expression[i]=='}'){ 
			if(openBrackets_curly.empty()) return false;
			if((!openBrackets_block.empty()) && openBrackets_block.back()>openBrackets_curly.back()) return false;
			if((!openBrackets_round.empty()) && openBrackets_round.back()>openBrackets_curly.back()) return false;
			pairedBrackets[openBrackets_curly.back()] = i; 
			pairedBrackets[i] = openBrackets_curly.back();
			openBrackets_curly.pop_back();
		}
	}
	if(!(openBrackets_curly.empty() && openBrackets_block.empty() && openBrackets_round.empty())) return false;
	return true;
}




bool MathExpression::findRoundBracketPairs(const string &expression, std::vector<long> &pairedBrackets){
	pairedBrackets.resize(expression.length(),-1);
	std::vector<long> openBrackets_round;
	for(long i=0; i<expression.length(); ++i){
		if(expression[i]=='('){ openBrackets_round.push_back(i); }
		else if(expression[i]==')'){ 
			if(openBrackets_round.empty()) return false;
			pairedBrackets[openBrackets_round.back()] = i; 
			pairedBrackets[i] = openBrackets_round.back();
			openBrackets_round.pop_back();
		}
	}
	if(!openBrackets_round.empty()) return false;
	return true;
}





bool MathExpression::parseAndEvaluate(	const string 	&expression, 		
										bool			allowRandoms,		
										bool			onlyRoundBrackets,
										string 			&errorMessage,
										double			&value){
	if(STMath_isReal(expression)){
		// expression is just a number
		value = STMath_string2Double(expression);
		return true;
	}else{
		// seems like a more complicated expression, so parse
		MathExpression ME;
		if(!ME.parse(expression, errorMessage, allowRandoms, onlyRoundBrackets, 0)) return false;
		value = ME();
		return true;
	}
}




bool MathExpression::parseAndEvaluate(	const string 	&expression, 		
										bool			allowRandoms,		
										bool			onlyRoundBrackets,
										string 			&errorMessage,
										double			&value,
										long			numberOfVariables,
										...){
	if(STMath_isReal(expression)){
		// expression is just a number
		value = STMath_string2Double(expression);
		return true;
	
	}else{
		// seems like a more complicated expression, so parse
		map<string, long> variableName2ID;
		std::vector<double> values(numberOfVariables);
		char *name_c;
		va_list vars;
		va_start(vars, numberOfVariables);
		for(long i=0; i<numberOfVariables; ++i){
			name_c = va_arg(vars, char*);
			values[i] = va_arg(vars, double);
			if(variableName2ID.find(name_c)!=variableName2ID.end()){ errorMessage = "Duplicate variable name '" + string(name_c) + "'"; return false; }
			variableName2ID[name_c] = i;
		}
		va_end(vars);
	
		// parse
		MathExpression ME;
		if(!ME.parse(expression, errorMessage, allowRandoms, onlyRoundBrackets, variableName2ID)){ errorMessage = "Could not parse expression: "+errorMessage; return false; }
		value = ME(values);
		return true;
	}
}





bool MathExpression::parse(	const string 			&expression, 
							string 					&errorMessage, 
							bool					allowRandoms,
							bool					onlyRoundBrackets,
							const map<string,long> 	&variableName2ID){
	return parse(expression, errorMessage, allowRandoms, onlyRoundBrackets, variableName2ID, map<string,long>());
}




bool MathExpression::parse(	const string 			&expression, 
							string 					&errorMessage,
							bool					allowRandoms,
							bool					onlyRoundBrackets,
							const map<string,long> 	&variableName2ID, 
							const map<string,long> 	&moreVariableName2ID){
	clear();
	originalExpression 	= expression;
	hasRandoms			= false;
	
	// figure out bracket pairs
	std::vector<long> pairedBrackets;
	if(onlyRoundBrackets){
		if(!findRoundBracketPairs(expression, pairedBrackets)){
			errorMessage = "Bracket mismatch"; return false;
		}
	}else{
		if(!findBracketPairs(expression, pairedBrackets)){
			errorMessage = "Bracket mismatch"; return false;
		}	
	}
	
	
	// prepare functional stack
	stackArguments.assign(1,std::vector<long>());
	stackValues.assign(1,0);
	stackFunction.assign(1,FunctionTypeUnknown);
			
	errorMessage = parseBlock(expression, 0, expression.length()-1, pairedBrackets, variableName2ID, moreVariableName2ID, 0, allowRandoms);
	if(!(OK = (errorMessage == ""))) return false;
	
	// if expression can be evaluated right away (i.e. is not random and does not depend on any variable), then optimize call stack by evaluating it right away
	if((!allowRandoms) && variableNames.empty()){
		for(long i=stackValues.size()-1; i>=0; --i){
			evaluateStackEntry(i);
		}
		setToConstant(stackValues[0]);
	}
	return true;
}



bool MathExpression::parse(	const string 	&expression, 
							string 			&errorMessage, 
							bool			allowRandoms,
							bool			onlyRoundBrackets,
							long 			numberOfVariables, 
							...){
	map<string, long> variableName2ID;
	
	// get variable names
	char *name_c;
	va_list vars;
	va_start(vars, numberOfVariables);
	for(long i=0; i<numberOfVariables; ++i){
		name_c = va_arg(vars, char*);
		if(variableName2ID.find(name_c)!=variableName2ID.end()){ errorMessage = "Duplicate variable name '" + string(name_c) + "'"; return false; }
		variableName2ID[name_c] = i;
	}
	va_end(vars);
	
	// parse
	return parse(expression, errorMessage, allowRandoms, onlyRoundBrackets, variableName2ID, map<string,long>());
}



bool MathExpression::parse(	const string 			&expression, 
							string 					&errorMessage, 
							bool					allowRandoms,
							bool					onlyRoundBrackets,
							long 					numberOfPrimaryVariables, 
							const map<string,long> 	&primaryVariableName2ID,
							long 					numberOfAdditionalVariables, 
							...){
				
	map<string, long> moreVariableName2ID;
	
	// get additional variable names
	char *name_c;
	va_list vars;
	va_start(vars, numberOfAdditionalVariables);
	for(long i=0; i<numberOfAdditionalVariables; ++i){
		name_c = va_arg(vars, char*);
		if((primaryVariableName2ID.find(name_c)!=primaryVariableName2ID.end()) || (moreVariableName2ID.find(name_c)!=moreVariableName2ID.end())){ errorMessage = "Duplicate variable name '" + string(name_c) + "'"; return false; }
		moreVariableName2ID[name_c] = numberOfPrimaryVariables + i;
	}
	va_end(vars);
	
	// parse
	return parse(expression, errorMessage, allowRandoms, onlyRoundBrackets, primaryVariableName2ID, moreVariableName2ID);
}





string MathExpression::functionType2str(FunctionType f){
	switch(f){
	case FunctionTypeRNormal: return "rnormal";
	case FunctionTypeRPoisson: return "rpoisson";
	case FunctionTypeRBernoulli: return "rbernoulli";
	case FunctionTypeRBinomial: return "rbinomial";
	case FunctionTypeRUniform: return "runiform";
	case FunctionTypeRLogUniform: return "rloguniform";
	case FunctionTypeRCauchy: return "rcauchy";
	case FunctionTypeRChiSquared: return "rchisquared";
	case FunctionTypeRTriangular: return "rtriangular";		
	case FunctionTypePlus: return "+";		
	case FunctionTypeMinus: return "-";	
	case FunctionTypeMultiply:	return "*";
	case FunctionTypeDivide: return "/";
	case FunctionTypeModulo: return "%";
	case FunctionTypeMin: return "min";
	case FunctionTypeMax: return "max";
	case FunctionTypeCos: return "cos";
	case FunctionTypeSin:	return "sin";
	case FunctionTypeTan:	return "tan";
	case FunctionTypeCot:	return "cot";
	case FunctionTypeAcos: return "acos";
	case FunctionTypeAsin: return "asin";
	case FunctionTypeAtan: return "atan";
	case FunctionTypeAtan2: return "atan2";
	case FunctionTypeAcot: return "acot";
	case FunctionTypeCosh: return "cosh";
	case FunctionTypeSinh: return "sinh";
	case FunctionTypeTanh: return "tanh";
	case FunctionTypeCoth: return "coth";
	case FunctionTypeExp: return "exp";
	case FunctionTypeLog: return "log";
	case FunctionTypeLog10: return "log10";
	case FunctionTypePow: return "^";
	case FunctionTypeSqrt: return "sqrt";
	case FunctionTypeCeil: return "ceil";
	case FunctionTypeFloor: return "floor";
	case FunctionTypeAbs: return "abs";
	case FunctionTypeHeaviside: return "heaviside";
	case FunctionTypePulse: return "pulse";
	case FunctionTypeComb: return "comb";
	case FunctionTypeEscapeNAN: return "escapeNAN";
	case FunctionTypeEscapeNAN2: return "escapeNAN2";
	case FunctionTypeEscapeInf: return "escapeInf";
	case FunctionTypeEscapeInf2: return "escapeInf2";
	case FunctionTypePiecewise2: return "piecewise2";
	case FunctionTypePiecewise3: return "piecewise3";
	case FunctionTypeNegate: return "negate";
	case FunctionTypeEvaluateVariable: return "variable";
	case FunctionTypeNumericalConstant: return "constant";
	case FunctionTypeUnknown: return "unknown";
	default: throw std::runtime_error("MathExpression: Unknown function type. Maybe a bug?");
	}
}




string MathExpression::functionType2description(FunctionType f){
	switch(f){
	case FunctionTypeRNormal: return "normal distribution";
	case FunctionTypeRPoisson: return "Poisson distribution";
	case FunctionTypeRBernoulli: return "Bernoulli distribution";
	case FunctionTypeRBinomial: return "binomial distribution";
	case FunctionTypeRUniform: return "uniform distribution";
	case FunctionTypeRLogUniform: return "log-uniform distribution";
	case FunctionTypeRCauchy: return "Cauchy distribution";
	case FunctionTypeRChiSquared: return "chi-squared distribution";
	case FunctionTypeRTriangular: return "triangular distribution";		
	case FunctionTypePlus: return "addition";		
	case FunctionTypeMinus: return "subtraction";	
	case FunctionTypeMultiply:	return "multiplication";
	case FunctionTypeDivide: return "division";
	case FunctionTypeModulo: return "modulo";
	case FunctionTypePow: return "exponentiation";
	case FunctionTypeMin: return "minimum";
	case FunctionTypeMax: return "maximum";
	case FunctionTypeCos: return "cosine";
	case FunctionTypeSin:	return "sine";
	case FunctionTypeTan:	return "tangent";
	case FunctionTypeCot:	return "cotangent";
	case FunctionTypeAcos: return "arccosinus";
	case FunctionTypeAsin: return "arcsinus";
	case FunctionTypeAtan: return "arctangent (slope-based)";
	case FunctionTypeAtan2: return "arctangent (coordinates-based)";
	case FunctionTypeAcot: return "arccotangent";
	case FunctionTypeCosh: return "hyberbolic cosine";
	case FunctionTypeSinh: return "hyberbolic sine";
	case FunctionTypeTanh: return "hyberbolic tangent";
	case FunctionTypeCoth: return "hyberbolic cotangent";
	case FunctionTypeExp: return "exponential";
	case FunctionTypeLog: return "natural logarithm";
	case FunctionTypeLog10: return "decadic logarithm";
	case FunctionTypeSqrt: return "square root";
	case FunctionTypeCeil: return "ceiling (next-highest integer)";
	case FunctionTypeFloor: return "floor (next-lowest integer)";
	case FunctionTypeAbs: return "absolute value";
	case FunctionTypeHeaviside: return "Heaviside step function";
	case FunctionTypePulse: return "rectangular pulse function";
	case FunctionTypeComb: return "periodic rectangular pulse function";
	case FunctionTypeEscapeNAN: return "NAN-conditional (1 escape choice)";
	case FunctionTypeEscapeNAN2: return "NAN-conditional (2 choices)";
	case FunctionTypeEscapeInf: return "INF-conditional (1 escape choice)";
	case FunctionTypeEscapeInf2: return "INF-conditional (2 choices)";
	case FunctionTypePiecewise2: return "piecewise function (2 choices)";
	case FunctionTypePiecewise3: return "piecewise function (3 choices)";
	case FunctionTypeNegate: return "negation";
	case FunctionTypeEvaluateVariable: return "variable";
	case FunctionTypeNumericalConstant: return "constant";
	case FunctionTypeUnknown: return "unknown";
	default: throw std::runtime_error("MathExpression: Unknown function type. Maybe a bug?");
	}
}




string MathExpression::functionType2genericExample(FunctionType f){
	switch(f){
	case FunctionTypeRNormal: return functionType2str(f)+"(mean,standard_deviation)";
	case FunctionTypeRPoisson: return functionType2str(f)+"(mean)";
	case FunctionTypeRBernoulli: return functionType2str(f)+"(p)";
	case FunctionTypeRBinomial: return functionType2str(f)+"(n,p)";
	case FunctionTypeRUniform: return functionType2str(f)+"(min,max)";
	case FunctionTypeRLogUniform: return functionType2str(f)+"(min,max)";
	case FunctionTypeRCauchy: return functionType2str(f)+"(median,scale)";
	case FunctionTypeRChiSquared: return functionType2str(f)+"(n)";
	case FunctionTypeRTriangular: return functionType2str(f)+"(mode,min,max)";		
	case FunctionTypePlus: return "a "+functionType2str(f)+" b";	
	case FunctionTypeMinus: return "a "+functionType2str(f)+" b";
	case FunctionTypeMultiply:	return "a "+functionType2str(f)+" b";
	case FunctionTypeDivide: return "a "+functionType2str(f)+" b";
	case FunctionTypeModulo: return "divident "+functionType2str(f)+" divisor";
	case FunctionTypePow: return "base "+functionType2str(f)+" exponent";
	case FunctionTypeMin: return functionType2str(f)+"(a,b)";
	case FunctionTypeMax: return functionType2str(f)+"(a,b)";
	case FunctionTypeCos: return functionType2str(f)+"(angle)";
	case FunctionTypeSin:	return functionType2str(f)+"(angle)";
	case FunctionTypeTan:	return functionType2str(f)+"(angle)";
	case FunctionTypeCot:	return functionType2str(f)+"(angle)";
	case FunctionTypeAcos: return functionType2str(f)+"(x)";
	case FunctionTypeAsin: return functionType2str(f)+"(x)";
	case FunctionTypeAtan: return functionType2str(f)+"(slope)";
	case FunctionTypeAtan2: return functionType2str(f)+"(y,x)";
	case FunctionTypeAcot: return functionType2str(f)+"(slope)";
	case FunctionTypeCosh: return functionType2str(f)+"(x)";
	case FunctionTypeSinh: return functionType2str(f)+"(x)";
	case FunctionTypeTanh: return functionType2str(f)+"(x)";
	case FunctionTypeCoth: return functionType2str(f)+"(x)";
	case FunctionTypeExp: return functionType2str(f)+"(x)";
	case FunctionTypeLog: return functionType2str(f)+"(x)";
	case FunctionTypeLog10: return functionType2str(f)+"(x)";
	case FunctionTypeSqrt: return functionType2str(f)+"(x)";
	case FunctionTypeCeil: return functionType2str(f)+"(x)";
	case FunctionTypeFloor: return functionType2str(f)+"(x)";
	case FunctionTypeAbs: return functionType2str(f)+"(x)";
	case FunctionTypeHeaviside: return functionType2str(f)+"(x)";
	case FunctionTypePulse: return functionType2str(f)+"(time,start,duration)";
	case FunctionTypeComb: return functionType2str(f)+"(time,start,duration,period)";
	case FunctionTypeNegate: return "-x";
	case FunctionTypeEscapeNAN: return functionType2str(f)+"(x,alternative_if_x_NAN)";
	case FunctionTypeEscapeNAN2: return functionType2str(f)+"(x,y_if_x_not_NAN,y_if_x_NAN)";
	case FunctionTypeEscapeInf: return functionType2str(f)+"(x,alternatie_if_x_INF)";
	case FunctionTypeEscapeInf2: return functionType2str(f)+"(x,y_if_x_not_INF,y_if_x_INF)";
	case FunctionTypePiecewise2: return functionType2str(f)+"(x,threshold,y_left,y_right)";
	case FunctionTypePiecewise3: return functionType2str(f)+"(x,threshold1,threshold2,y_left,y_middle,y_right)";
	case FunctionTypeEvaluateVariable: return "x";
	case FunctionTypeNumericalConstant: return "0";
	case FunctionTypeUnknown: return "unknown";
	default: throw std::runtime_error("MathExpression: Unknown function type. Maybe a bug?");
	}
}




long MathExpression::functionType2dim(FunctionType f){
	switch(f){
	case FunctionTypePiecewise2: return 4;
	case FunctionTypePiecewise3: return 6;
	case FunctionTypeComb: return 4;
	case FunctionTypeRTriangular: return 3;		
	case FunctionTypeRNormal: return 2;
	case FunctionTypeRBinomial: return 2;
	case FunctionTypeRUniform: return 2;
	case FunctionTypeRLogUniform: return 2;
	case FunctionTypeRCauchy: return 2;
	case FunctionTypePlus:	 return 2;	
	case FunctionTypeMinus: return 2;
	case FunctionTypeMultiply: return 2;
	case FunctionTypeDivide: return 2;
	case FunctionTypePow: return 2;
	case FunctionTypeAtan2: return 2;
	case FunctionTypeMin: return 2;
	case FunctionTypeMax: return 2;
	case FunctionTypeModulo: return 2;
	case FunctionTypePulse: return 3;
	case FunctionTypeRChiSquared: return 1;
	case FunctionTypeRBernoulli: return 1;
	case FunctionTypeRPoisson: return 1;
	case FunctionTypeCos: return 1;
	case FunctionTypeSin: return 1;
	case FunctionTypeTan: return 1;
	case FunctionTypeCot: return 1;
	case FunctionTypeAcos: return 1;
	case FunctionTypeAsin: return 1;
	case FunctionTypeAtan: return 1;
	case FunctionTypeAcot: return 1;
	case FunctionTypeCosh: return 1;
	case FunctionTypeSinh: return 1;
	case FunctionTypeTanh: return 1;
	case FunctionTypeCoth: return 1;
	case FunctionTypeExp: return 1;
	case FunctionTypeLog: return 1;
	case FunctionTypeLog10: return 1;
	case FunctionTypeSqrt: return 1;
	case FunctionTypeCeil: return 1;
	case FunctionTypeFloor: return 1;
	case FunctionTypeAbs: return 1;
	case FunctionTypeHeaviside: return 1;
	case FunctionTypeNegate: return 1;
	case FunctionTypeEscapeNAN: return 2;
	case FunctionTypeEscapeNAN2: return 3;
	case FunctionTypeEscapeInf: return 2;
	case FunctionTypeEscapeInf2: return 3;
	case FunctionTypeEvaluateVariable: return 0;
	case FunctionTypeNumericalConstant: return 0;
	default: throw std::runtime_error("MathExpression: Unknown function type");
	}
}



bool MathExpression::functionTypeRandom(FunctionType f){
	switch(f){
	case FunctionTypeRNormal:
	case FunctionTypeRPoisson:
	case FunctionTypeRBernoulli:
	case FunctionTypeRBinomial:
	case FunctionTypeRUniform:
	case FunctionTypeRLogUniform:
	case FunctionTypeRCauchy:
	case FunctionTypeRChiSquared:
	case FunctionTypeRTriangular: return true;		
	default: return false;
	}
}


bool MathExpression::functionTypeOperator(FunctionType f){
	switch(f){
	case FunctionTypePlus:
	case FunctionTypeMinus:
	case FunctionTypeMultiply:
	case FunctionTypeDivide:
	case FunctionTypeModulo:
	case FunctionTypePow: return true;
	default: return false;
	}
}



void MathExpression::printDebug(ostream &stream) const{
	stream << "Math expression is " << (OK ? "OK" : "not OK") << "\n  Stack contains " << variableValues.size() << " variables and " << stackValues.size() << " operations\n";
	stream << "  Variables:\n";
	for(long i=0; i<variableValues.size(); ++i){
		stream << "  " << i << ":" << variableNames[i] << " (ID " << variableIDs[i] << ") = " << stackValues[i] << "\n";
	}
	stream << "  Operations:\n";
	for(long i=0; i<stackValues.size(); ++i){
		if(stackFunction[i]==FunctionTypeEvaluateVariable){ stream << "  " << i << ":" << (stackArguments[i][0]>=0 ? variableNames[stackArguments[i][0]] : " NA") << "\n"; }
		else{ 
			stream << "  " << i << ":" << functionType2str(stackFunction[i]) << " ("; 
			printTuple(stream, stackArguments[i]);
			stream << ") = " << stackValues[i] << "\n"; 
		}
	}
}






double MathExpression::operator()(const std::vector<double> &x) const{
	if(stackValues.empty()) return 0;
	if(!OK) return 0;
	if((stackFunction.size()==1) && (stackFunction[0]==FunctionTypeNumericalConstant)) return stackValues[0];
	const long NV = variableNames.size();
	long i;
	
	// extract relevant variable values
	for(i=0; i<NV; ++i){
		if(variableIDs[i]>=x.size()){ variableValues[i]=0; }
		else{ variableValues[i]=x[variableIDs[i]]; }
	}

	// evaluate the evaluation stack
	for(i=stackValues.size()-1; i>=0; --i){
		evaluateStackEntry(i);
	}

	// return root entry (last calculated)
	return stackValues[0];

}




double MathExpression::evaluateAt(long numberOfVariables, ...) const{
	if(stackValues.empty()) return 0;
	const long NV = variableNames.size();
	if(!OK) return 0;
	if((stackFunction.size()==1) && (stackFunction[0]==FunctionTypeNumericalConstant)) return stackValues[0];
	double x;
	long i;

	// extract relevant variable values
	va_list vars;
	va_start(vars, numberOfVariables);
	map<long,long>::const_iterator id2v;
	for(i=0; i<numberOfVariables; ++i){
		x = va_arg(vars, double);
		if((id2v=ID2variable.find(i))!=ID2variable.end()){
			variableValues[id2v->second] = x;
		}
	}
	va_end(vars);
	
	// set missing variables to zero
	for(i=0; i<NV; ++i){
		if(variableIDs[i]>=numberOfVariables){ variableValues[i]=0; }
	}
	
	// evaluate the evaluation stack
	for(i=stackValues.size()-1; i>=0; --i){
		evaluateStackEntry(i);
	}
	
	// return root entry (last calculated)
	return stackValues[0];
}



double MathExpression::operator()(const std::vector<double> &x, long numberOfAdditionalVariables, ...) const{
	if(stackValues.empty()) return 0;
	if(!OK) return 0;
	if((stackFunction.size()==1) && (stackFunction[0]==FunctionTypeNumericalConstant)) return stackValues[0];
	const long NV = variableNames.size();
	long i;
	double ax;
	
	// extract relevant variable values
	for(i=0; i<NV; ++i){
		if(variableIDs[i]>=x.size()+numberOfAdditionalVariables){ variableValues[i]=0; }
		else if(variableIDs[i]<x.size()){ variableValues[i]=x[variableIDs[i]]; }
	}
	va_list vars;
	va_start(vars, numberOfAdditionalVariables);
	map<long,long>::const_iterator id2v;
	for(i=0; i<numberOfAdditionalVariables; ++i){
		ax = va_arg(vars, double);
		if((id2v=ID2variable.find(x.size()+i))!=ID2variable.end()){
			variableValues[id2v->second] = ax;
		}
	}
	va_end(vars);

	// evaluate the evaluation stack
	for(i=stackValues.size()-1; i>=0; --i){
		evaluateStackEntry(i);
	}

	// return root entry (last calculated)
	return stackValues[0];
}



double MathExpression::operator()(const std::vector<double> &x1, const std::vector<double> &x2, long numberOfAdditionalVariables, ...) const{
	if(stackValues.empty()) return 0;
	if(!OK) return 0;
	if((stackFunction.size()==1) && (stackFunction[0]==FunctionTypeNumericalConstant)) return stackValues[0];
	const long NV = variableNames.size();
	const long N1 = x1.size();
	const long N2 = x2.size();
	const long N3 = numberOfAdditionalVariables;
	long i;
	double ax;
	
	// extract relevant variable values
	for(i=0; i<NV; ++i){
		if(variableIDs[i]>=N1+N2+N3){ variableValues[i]=0; }
		else if(variableIDs[i]<N1){ variableValues[i]=x1[variableIDs[i]]; }
		else if(variableIDs[i]<N1+N2){ variableValues[i]=x2[variableIDs[i]-N1]; }
	}
	va_list vars;
	va_start(vars, numberOfAdditionalVariables);
	map<long,long>::const_iterator id2v;
	for(i=0; i<N3; ++i){
		ax = va_arg(vars, double);
		if((id2v=ID2variable.find(N1+N2+i))!=ID2variable.end()){
			variableValues[id2v->second] = ax;
		}
	}
	va_end(vars);

	// evaluate the evaluation stack
	for(i=stackValues.size()-1; i>=0; --i){
		evaluateStackEntry(i);
	}

	// return root entry (last calculated)
	return stackValues[0];
}



double MathExpression::operator()(const std::vector<double> &x1, const std::vector<double> &x2, const std::vector<double> &x3, long numberOfAdditionalVariables, ...) const{
	if(stackValues.empty()) return 0;
	if(!OK) return 0;
	if((stackFunction.size()==1) && (stackFunction[0]==FunctionTypeNumericalConstant)) return stackValues[0];
	const long NV = variableNames.size();
	const long N1 = x1.size();
	const long N2 = x2.size();
	const long N3 = x3.size();
	const long N4 = numberOfAdditionalVariables;
	long i;
	double ax;
	
	// extract relevant variable values
	for(i=0; i<NV; ++i){
		if(variableIDs[i]>=N1+N2+N3+N4){ 	variableValues[i] = 0; }
		else if(variableIDs[i]<N1){ 		variableValues[i] = x1[variableIDs[i]]; }
		else if(variableIDs[i]<N1+N2){ 		variableValues[i] = x2[variableIDs[i]-N1]; }
		else if(variableIDs[i]<N1+N2+N3){ 	variableValues[i] = x3[variableIDs[i]-N1-N2]; }
	}
	va_list vars;
	va_start(vars, numberOfAdditionalVariables);
	map<long,long>::const_iterator id2v;
	for(i=0; i<N4; ++i){
		ax = va_arg(vars, double);
		if((id2v=ID2variable.find(N1+N2+N3+i))!=ID2variable.end()){
			variableValues[id2v->second] = ax;
		}
	}
	va_end(vars);

	// evaluate the evaluation stack
	for(i=stackValues.size()-1; i>=0; --i){
		evaluateStackEntry(i);
	}

	// return root entry (last calculated)
	return stackValues[0];
}



double MathExpression::operator()(long numberOfVariableVectors, long numberOfSingleVariables, ...) const{
	if(stackValues.empty()) return 0;
	if(!OK) return 0;
	if((stackFunction.size()==1) && (stackFunction[0]==FunctionTypeNumericalConstant)) return stackValues[0];
	const long NV = variableNames.size();
	long cumulative;
	std::vector<std::vector<double>* > variableVectors(numberOfVariableVectors);
	
	// extract relevant variable values
	long Nprovided = numberOfSingleVariables;
	va_list vars;
	va_start(vars, numberOfSingleVariables);
	for(long i=0; i<numberOfVariableVectors; ++i){
		variableVectors[i] = va_arg(vars, std::vector<double>*);
		Nprovided += variableVectors[i]->size();
	}
	for(long i=0; i<NV; ++i){
		if(variableIDs[i]>=Nprovided){ variableValues[i] = 0; continue; }
		for(long j=cumulative=0; j<numberOfVariableVectors; ++j){
			cumulative += variableVectors[j]->size();
			if(variableIDs[i]<cumulative){ variableValues[i]=(*variableVectors[j])[variableIDs[i]-cumulative+variableVectors[j]->size()]; break; }
		}
	}
	map<long,long>::const_iterator id2v;
	double ax;
	for(long j=Nprovided-numberOfSingleVariables; j<Nprovided; ++j){
		ax = va_arg(vars, double);
		if((id2v=ID2variable.find(j))!=ID2variable.end()){
			variableValues[id2v->second] = ax;
		}
	}
	va_end(vars);
	
	// evaluate the evaluation stack
	for(long i=stackValues.size()-1; i>=0; --i){
		evaluateStackEntry(i);
	}

	// return root entry (last calculated)
	return stackValues[0];
}




bool MathExpression::dependsOnVariable(const string &name) const{
	for(long i=0; i<variableNames.size(); ++i){
		if(variableNames[i]==name) return true;
	}
	return false;
}

bool MathExpression::dependsOnVariableWithPrefix(const string &namePrefix) const{
	for(long i=0; i<variableNames.size(); ++i){
		if(variableNames[i].substr(0, namePrefix.length()) == namePrefix) return true;
	}
	return false;	
}


void MathExpression::setToConstant(double C){
	variableNames.clear();
	variableValues.clear();
	ID2variable.clear();
	variableIDs.clear();
	stackValues.assign(1,C);
	stackFunction.assign(1,FunctionTypeNumericalConstant);
	stackArguments.assign(1,std::vector<long>());
	originalExpression = STMath_makeString(C);
	OK = true;
	hasRandoms = false;
}




template<class TYPE>
void MathExpression::printTuple(ostream &sout, const std::vector<TYPE> &components) const{
	for(long i=0; i<components.size(); ++i){
		sout << (i==0 ? "" : ", ") << components[i];
	}
}


void MathExpression::getListOfAvailableFunctions(	std::vector<string> &names, 
													std::vector<long> 	&numberOfArguments,
													std::vector<string> &descriptions,
													std::vector<string>	&genericExamples,
													bool 				includeRandom, 
													bool 				includeDeterministic, 
													bool 				includeOperators){
	names.clear(); numberOfArguments.clear();
	for(FunctionType f=FunctionType(FunctionTypeStartOfFunctionList+1); f<FunctionTypeEndOfFunctionList; f = FunctionType(f+1)){
		if(functionTypeRandom(f)){
			if(!includeRandom) continue;
		}else if(functionTypeOperator(f)){
			if(!includeOperators) continue;
		}else{
			if(!includeDeterministic) continue;
		}
		names.push_back(functionType2str(f));
		numberOfArguments.push_back(functionType2dim(f));
		descriptions.push_back(functionType2description(f));
		genericExamples.push_back(functionType2genericExample(f));
	}
}


bool MathExpression::isValidName(const string &s, string &errorMessage){
	for(long i=0; i<s.size(); ++i){
		if(!isValidNameCharacter(s[i])){
			errorMessage = "Invalid character '" + s.substr(i,1) + "'";
			return false;
		}
	}
	return true;
}


string MathExpression::lowercase(string s){
	for(long n=0; n<s.length(); ++n){ s[n] = tolower(s[n]); }
	return s;
}




double MathExpression::aux_pulse(double time, double start, double duration){
	return (time<start || time>start+duration ? 0 : 1);
}


double MathExpression::aux_comb(double time, double start, double duration, double periodicity){
	return aux_pulse(fmod(time-start, periodicity), 0.0, duration);
}



double MathExpression::random_uniformWithinInclusiveRight(double minimum, double maximum){
	return minimum + (maximum - minimum)*R::runif(0.0+1e30,1.0);
}

double MathExpression::random_uniformWithinInclusiveLeft(double minimum, double maximum){
	return minimum + (maximum - minimum)*R::runif(0.0,1.0-1e-30);
}

double MathExpression::random_logUniformWithinInclusive(double minimum, double maximum){
	return exp(R::runif(log(minimum), log(maximum)));
}

double MathExpression::random_standardNormal(){
	return sqrt(-2.0*log(random_uniformWithinInclusiveRight(0, 1)))*cos(2.0*M_PI*random_uniformWithinInclusiveRight(0,1));
}

double MathExpression::random_normal(double mean, double standardDeviation){
	return random_standardNormal()*standardDeviation + mean;
}

bool MathExpression::random_bernoulli(double p){ 
	return (R::runif(0.0,1.0)<p);
}

long MathExpression::random_binomial(long N, double p){
	long count = 0;
	for(long n=0; n<N; ++n){
		count += (random_bernoulli(p) ? 1 : 0);
	}
	return count;
}

long MathExpression::random_Poisson(double mean){
	long k=0;
	double L = exp(-mean);
	if(mean<50){
		double p = R::runif(0.0,1.0);
		double cum=L;
		while(cum<p){
			++k;
			L *= mean/k;
			cum += L;
		}
		return k;
	}else{
		//approximate by normal distribution
		return max(0l, (long)(random_standardNormal()*sqrt(mean)+mean));
	}
}

inline double MathExpression::random_Cauchy(double mean, double gamma){
	return mean + gamma * tan(M_PI * (R::runif(0.0,1.0) - 0.5));
}

double MathExpression::random_chiSquare(long degrees){
	double x=0, dummy; long i;
	for(i=0; i<degrees; ++i){
		dummy = random_standardNormal();
		x += dummy*dummy;
	}
	return x;
}

//returns triangularly-distributed number. 
double MathExpression::random_triangular(double mode, double triangleMin, double triangleMax){
	double x = R::runif(0.0,1.0);
	mode = max(triangleMin, min(triangleMax, mode));
	if(x*(triangleMax - triangleMin) < (mode - triangleMin)){
		return triangleMin + sqrt(x * (triangleMax - triangleMin)*(mode - triangleMin));
	}else{
		return triangleMax - sqrt((1-x)*(triangleMax - triangleMin)*(triangleMax - mode));
	}
}





// evaluate a univariate mathematical expression (function of some predictor variable X, C-syntax) for multiple input X-values
// [[Rcpp::export]]
Rcpp::List evaluate_univariate_expression_CPP(	const std::string 			&expression,	// valid mathematical expression, depending on a single variable
												const std::string 			&Xname,			// name of the X-variable, i.e. as it appears in the expression. For example "$x" or "$x$".
												const std::vector<double>	&X){			// 1D array of X-values, at which to evaluate expression
	// parse expression into a mathematical object
	MathExpression parser;
	string error;
	parser.parse(expression, error, true, false, 1, Xname.c_str());
	if(error!="") return Rcpp::List::create(Rcpp::Named("success") = false, Rcpp::Named("error") = error);

	// evaluate expression for each X
	std::vector<double> Y(X.size());
	for(long i=0; i<X.size(); ++i){
		Y[i] = parser(X[i]);
	}

	return Rcpp::List::create(Rcpp::Named("success")  = true, Rcpp::Named("Y") = Y);
}



#endif

