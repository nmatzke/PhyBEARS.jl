// 	STMathExpression
//	Class for evaluating mathematical expressions on runtime
//	Includes support for random numbers and symbolic variables in math expressions
//	
//	Stilianos Louca
//	August 31, 2015


#ifndef STMATHEXPRESSION_DECL
#define STMATHEXPRESSION_DECL

#include <new>
#include <Rcpp.h>
#include <vector>
#include <limits>
#include <stdexcept>
#include <string>
#include <map>

#ifndef INFTY_D
#define INFTY_D numeric_limits<double>::infinity()
#endif

using namespace std;

typedef enum {
	// formal placeholder for beginning of functions available in a MathExpression. This needs to be at the start of the function list
	FunctionTypeStartOfFunctionList,

	// random functions
	FunctionTypeRNormal,
	FunctionTypeRPoisson,
	FunctionTypeRBernoulli,
	FunctionTypeRBinomial,
	FunctionTypeRUniform,
	FunctionTypeRLogUniform,
	FunctionTypeRCauchy,
	FunctionTypeRChiSquared,
	FunctionTypeRTriangular,

	// binary operators
	FunctionTypePlus,
	FunctionTypeMinus,
	FunctionTypeMultiply,
	FunctionTypeDivide,
	FunctionTypeModulo,
	FunctionTypePow,

	// deterministic functions
	FunctionTypeMin,
	FunctionTypeMax,
	FunctionTypeCos,
	FunctionTypeSin,
	FunctionTypeTan,
	FunctionTypeCot,
	FunctionTypeAcos,
	FunctionTypeAsin,
	FunctionTypeAtan,
	FunctionTypeAtan2,
	FunctionTypeAcot,
	FunctionTypeCosh,
	FunctionTypeSinh,
	FunctionTypeTanh,
	FunctionTypeCoth,
	FunctionTypeExp,
	FunctionTypeLog,
	FunctionTypeLog10,
	FunctionTypeSqrt,
	FunctionTypeCeil,
	FunctionTypeFloor,
	FunctionTypeAbs,
	FunctionTypeHeaviside,
	FunctionTypePulse,
	FunctionTypeComb,
	FunctionTypeEscapeNAN,				// escapeNAN(x,y) evaluates to y if x==NAN, and to x otherwise
	FunctionTypeEscapeNAN2,				// escapeNAN2(x,y1,y2) evaluates to y1 if x!=NAN, and to y2 if x==NAN
	FunctionTypeEscapeInf,				// escapeInf(x,y) evaluates to y if x==+-INFTY_D, and to x otherwise
	FunctionTypeEscapeInf2,				// escapeInf2(x,y,z) evaluates to y if x==+INFTY_D, to z if x=-INFTY_D, and to x otherwise
	FunctionTypePiecewise2, 			// piecewise2(x,a,y1,y2) = (x<a ? y1 : y2)
	FunctionTypePiecewise3,				// piecewise3(x,a,b,y1,y2,y3) = (x<a ? y1 : (x<b ? y2 : y3))
	
	// formal placeholder for end of functions available in a MathExpression. This needs to be at the end of the function list
	FunctionTypeEndOfFunctionList,
	
	// internal (formal) functions, for the logistics of MathExpression
	FunctionTypeNegate,
	FunctionTypeEvaluateVariable,
	FunctionTypeNumericalConstant,
	FunctionTypeUnknown					// placeholder for undefined function
} FunctionType;








// structure for parsing and evaluating arbitrary mathematical expressions/functions of arbitrary number of variables
// parsing is computationally expensive, but evaluation is pretty fast (since call-stack is allready set up)
// function names (but not variables or constants) are case-insensitive
// defaults to zero
class MathExpression{
private:
	std::vector<std::vector<long> >	stackArguments;
	mutable std::vector<double>		stackValues;
	std::vector<FunctionType>		stackFunction;
	string	 						originalExpression;
	
	mutable std::vector<double>		variableValues;
	std::vector<string>				variableNames;
	std::vector<long>				variableIDs;
	std::map<long,long>				ID2variable;
	
	bool	OK;
	bool	hasRandoms;

	// some elementary math functions
	static double aux_pulse(double time, double start, double duration);
	static double aux_comb(double time, double start, double duration, double periodicity);
	
	// random number generators

	static double	random_uniformWithinInclusiveRight(double minimum, double maximum);
	static double 	random_uniformWithinInclusiveLeft(double minimum, double maximum);
	static double 	random_logUniformWithinInclusive(double minimum, double maximum);
	static double 	random_standardNormal();
	static double 	random_normal(double mean, double standardDeviation);
	static bool 	random_bernoulli(double p);
	static long 	random_binomial(long N, double p);
	static long 	random_Poisson(double mean);
	static double 	random_Cauchy(double mean, double gamma);
	static double 	random_chiSquare(long degrees);
	static double 	random_triangular(double mode, double triangleMin, double triangleMax);

	
	// parsing and evaluating
	void clear();
	void evaluateStackEntry(long i) const;
	
	// returns either "" on success or an error message
	string parseBlock(	const string 				&expression, 
						long 						start, 
						long 						end, 
						const std::vector<long> 	&pairedBrackets, 
						const std::map<string,long> &variableName2ID, 
						const std::map<string,long> &moreVariableName2ID, 
						long 						stackIndex, 
						bool 						allowRandoms);
	
	static FunctionType str2FunctionType(const string &s);
	static string functionType2str(FunctionType f);
	static string functionType2description(FunctionType f);
	static string functionType2genericExample(FunctionType f);
	static long functionType2dim(FunctionType f); // returns dimension of input tuple required by the function, i.e. number of arguments
	static bool functionTypeRandom(FunctionType f);
	static bool functionTypeOperator(FunctionType f);
	static string lowercase(string s);
	static bool is_one_of(char c, const char *list);
	
	static long getLeftMostOfPlusMinusOperatorChain(const string &expression, long start, long pos);
	static int 	getSignOfPlusMinusChain(const string &expression, long pos, long end, long &endOfChain);
	static bool findBracketPairs(const string &expression, std::vector<long> &pairedBrackets);
	static bool findRoundBracketPairs(const string &expression, std::vector<long> &pairedBrackets);
	static bool isWhiteSpace(const string &expression, long start, long end);
	static bool hasWhiteSpace(const string &expression, long start, long end);
	static long splitBinaryOperatorInverseOrder(const string &expression, long start, long end, const std::vector<long> &pairedBrackets);
	static long splitTuple(const string &expression, long start, long end, const std::vector<long> &pairedBrackets, std::vector<long> &commas);
	static long getBracketPrefix(const string &expression, long start, long br);
	bool parse(	const 					string &expression, 
				string 					&errorMessage, 
				bool					allowRandoms,
				bool					onlyRoundBrackets,
				const map<string,long> 	&variableName2ID, 
				const map<string,long> 	&moreVariableName2ID);
	
	template<class TYPE>
	void printTuple(ostream &sout, const std::vector<TYPE> &components) const;
	

public:
	MathExpression(){ clear(); }
	
	// parse mathematical expression
	// the ... arguments should be (const char *) types (variable names)
	// returns false on error
	// variable names must satisfy isValidName(..), and in particular only contain characters allowed by isValidNameCharacter(char c)
	bool parse(	const string 	&expression, 
				string 			&errorMessage, 
				bool			allowRandoms,	// if false, random functions are not allowed in expression
				bool			onlyRoundBrackets, // if false, all three bracket types "()", "{}", "[]" are allowed
				long 			numberOfVariables, 
				...);	// the remaining numberOfVariables variable names (const *char)
	
	// parse mathematical expression
	// variableName2ID defines the order and names of potential variables
	// variable names must satisfy isValidName(..), and in particular only contain characters allowed by isValidNameCharacter(char c)
	// returns false on error
	bool parse(	const string 			&expression, 
				string 					&errorMessage, 
				bool					allowRandoms,	// if false, random functions are not allowed in expression
				bool					onlyRoundBrackets, // if false, all three bracket types "()", "{}", "[]" are allowed
				const map<string,long> 	&variableName2ID);

	// similar to the previous parse, but with the possibility of defining further variable names
	// the remaining numberOfAdditionalVariables arguments should be of type (const char *)
	// variable names must satisfy isValidName(..), and in particular only contain characters allowed by isValidNameCharacter(char c)
	bool parse(	const string 			&expression, 
				string 					&errorMessage, 
				bool					allowRandoms,	// if false, random functions are not allowed in expression
				bool					onlyRoundBrackets, // if false, all three bracket types "()", "{}", "[]" are allowed
				long 					numberOfPrimaryVariables, 
				const map<string,long> 	&primaryVariableName2ID, // mapping from variable name to its ID (0,..,numberOfPrimaryVariables-1)
				long 					numberOfAdditionalVariables, 
				...);


	// parse and evaluate mathematical expression without variables
	static bool parseAndEvaluate(	const string 	&expression, 		// (INPUT) mathematical expression without any variables
									bool			allowRandoms,		// (INPUT) if false, random functions are not allowed in expression
									bool			onlyRoundBrackets, 	// (INPUT) if false, all three bracket types "()", "{}", "[]" are allowed
									string 			&errorMessage, 		// (OUTPUT) error message in case of an error
									double			&value);			// (OUTPUT) the value of the expression, in case of no error
	
	// parse and evaluate mathematical expression, with optional variables					
	// variable names must satisfy isValidName(..), and in particular only contain characters allowed by isValidNameCharacter(char c)
	static bool parseAndEvaluate(	const string 	&expression, 		// (INPUT) mathematical expression without any variables
									bool			allowRandoms,		// (INPUT) if false, random functions are not allowed in expression
									bool			onlyRoundBrackets, 	// (INPUT) if false, all three bracket types "()", "{}", "[]" are allowed
									string 			&errorMessage, 		// (OUTPUT) error message in case of an error
									double			&value,				// (OUTPUT) the value of the expression, in case of no error
									long 			numberOfVariables, 	// (INPUT) number of variables passed to the function. the remaining numberOfVariables arguments must be of alternating types (const *char) and double
									...);

	// evaluate parsed mathematical expression
	double operator()() const{ return evaluateAt(0); } // shortcut for constant expression evaluation (or setting all variables zero)
	double operator()(double x) const{ return evaluateAt(1,x); } // shortcut for univariate functions
	double operator()(const std::vector<double> &x) const; // values in x[] should correspond to original variables. If x is smaller than the original number of variables (as specified during parse), the rest is set to zero
	double operator()(const std::vector<double> &x, long numberOfAdditionalVariables, ...) const; // values in x[] should correspond to the first x.size() original variables. If the total number of passed variables is less than the original (as specified during parse), the rest is assumed zero
	double operator()(const std::vector<double> &x1, const std::vector<double> &x2, long numberOfAdditionalVariables, ...) const; // values in x1[] should correspond to the first x1.size() original variables, values in x2[] in the next x2.size() variables. "..." can contain any of the remaining variables. If the total number of passed variables is less than the original (as specified during parse), the rest is assumed zero
	double operator()(const std::vector<double> &x1, const std::vector<double> &x2, const std::vector<double> &x3, long numberOfAdditionalVariables, ...) const; // values in x1[] should correspond to the first x1.size() original variables, values in x2[] in the next x2.size() variables, values in x3[] in the next x3.size() variables. "..." can contain any of the remaining variables. If the total number of passed variables is less than the original (as specified during parse), the rest is assumed zero
	double operator()(long numberOfVariableVectors, long numberOfSingleVariables, ...) const; // "..." must contain numberOfVariableVectors vector<double>* pointers, followed by numberOfSingleVariables double. The passed vectors (VV1[], VV2[], ...) correspond to the first VV1.size + VV2.size + ... original variables. The remaining passed doubles (D1, D2, ...) correspond to the remaining original variables. If the total number of passed variables is less than the original (as specified during parse), the missing values are assumed zero
	double evaluateAt(long numberOfVariables, ...) const; // missing arguments should be double (variable values). If numberOfVariables is smaller than specified during parse, the remaining variables are assumed zero


	string getOriginalExpression() const{ return originalExpression; }
	std::vector<string> getVariableNames() const{ return variableNames; }
	bool dependsOnVariable(const string &name) const;
	bool dependsOnVariableWithPrefix(const string &namePrefix) const;
	bool isOK() const{ return OK; }
	bool isZero() const{ return (!OK) || stackValues.empty(); }
	void setToZero(){ clear(); }
	void setToConstant(double C);
	bool isConstant() const{ return (!(hasRandoms) && variableNames.empty()); }

	// prints cryptic debug message
	// for debugging only
	void printDebug(ostream &stream) const;
	
	inline static bool isValidNameCharacter(char c){ return (isalnum(c) || is_one_of(c, "_.@$!?")); }
	static bool isValidName(const string &s, string &errorMessage);

	static void getListOfAvailableFunctions(std::vector<string> 	&names, 
											std::vector<long> 	&numberOfArguments, 
											std::vector<string> &descriptions,
											std::vector<string>	&genericExamples,
											bool 				includeRandom, 
											bool 				includeDeterministic, 
											bool 				includeOperators);
};



// evaluate a univariate mathematical expression (function of some predictor variable X, C-syntax) for multiple input X-values
// [[Rcpp::export]]
Rcpp::List evaluate_univariate_expression_CPP(	const std::string 			&expression,	// valid mathematical expression, depending on a single variable
												const std::string 			&Xname,			// name of the X-variable, i.e. as it appears in the expression. For example "$x" or "$x$".
												const std::vector<double>	&X);			// 1D array of X-values, at which to evaluate expression








#endif

