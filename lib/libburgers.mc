/* [wxMaxima batch file version 1] [ DO NOT EDIT BY HAND! ]*/
/* [ Created with wxMaxima version 15.08.2 ] */

/* [wxMaxima: title   start ]
Librairie : Expectation function
   [wxMaxima: title   end   ] */

/* [wxMaxima: input   start ] */
/*  Utilitaire  */
range(imax):=makelist(i,i,1,imax)$
validationtest(test,funcname):=
/*  input
    -----
        test        :   list of boolean that represents result of elementary 
                            validation .
        funcname    :   name of function to validate
*/
block(
    [failure:false],
    for res in test do
    (
        if res=false then    
            failure:true
    ),        
    if failure then
        (
            print(test),
            throw("Error in validation of ** ",funcname," **")
        )    
    else
        print("Validation successfull ** ",funcname," **")
    )$
/* [wxMaxima: input   end   ] */

/* [wxMaxima: section start ]
Implémentation et validation de l'opérateur d'espérance
   [wxMaxima: section end   ] */

/* [wxMaxima: comment start ]
Dans cette section on définit et valide l'opérateur d'espérance qui intervient dans les 
calcul d'évolutio n des statistiques. 

To do: ne prend pas en compte les dérivées en temps.. d/dt
    --> Il faut traiter l'involution: E(E(c))=E(c)
   [wxMaxima: comment end   ] */

/* [wxMaxima: subsect start ]
Implementation de l'opérateur
   [wxMaxima: subsect end   ] */

/* [wxMaxima: comment start ]
L'opérateur E() repose sur les règles de calcules suivantes:
    
    (1) L'idempotence:  EoE=E
    (2) La linéarité:
        (a) E(X+Y)=E(X)+E(Y)
        (b) E(a*X)=a*E(X)       où 'a' est certaine.
    (3) La commutativité avec les opérateurs de différenciation d/dx, d/dt


L'implémentation de l'opérateur repose sur la distinction entre variable 
aléatoire et non-aléatoire. Ainsi on implémente dans un premier temps le
test si une entrée est considérée comme aléatoire ou non, puis la routine
E() est décrite.

    1) Test si une expression est aléatoire
    2) Reconnaissance automatique d'un motif: séparation a*X en partie aléatoire/non-aléa.
    3) Routine E()
   [wxMaxima: comment end   ] */

/* [wxMaxima: subsubsect start ]
Test de dépendance et test si une expression est aléatoire (israndom)
   [wxMaxima: subsubsect end   ] */

/* [wxMaxima: input   start ] */
/*

        ZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
                   Warning
        ZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
    
        Run the script until E() so the 
        routine is defined and can be used
        to tests within israndom.

*/


/**************************************  
Function:   **   israndom   **  
 **************************************/
/*  israndom(expr) : test if an expression expr is a function of ω in Ω the 
                probability space  
                a random variable is a function defined by
                defines(expr,ω)
                or a function f(ω)
*/
israndom(expr):=block(
                    /*
                    Input
                    -----
                    expr   :   expression
                    */
                        /*  Local variables         */
                        [output, a, opE],
                        output:false,
                        opE:op('E(a)),
                        /*  Core of the routine     */
                        if atom(expr) then
                            if expr=ω then
                                /*  case where expr=ω               */
                                output:true            
                            else
                                /*  Search of dependency            */
                                for fnc in dependencies do             
                                (                                    
                                    if op(fnc)=expr then
                                        output: not freeof(ω,fnc)
                                )
                        else
                            /*  When f is a function given by its arguments */
                            (
                                /*  E() is certain... so end */
                                if is(op(expr)#opE) then
                                    (
                                    /*   (i) test if 'ω' is an argument of f    */
                                    /*          case of ε(ω) for random var.    */
                                    output: not freeof(ω,expr),
                                    /*  (ii) test if 'ω' is implicitly an 
                                            argument of f                   */
                                    /*          case of f(ε) for random var.
                                                since equal to f[ε(ω)]      */
                                    if output=false then
                                        for arg in args(expr) do
                                            (
                                                if israndom(arg) then
                                                    output:true                                                                 
                                            ),
                                    /*  (iii) test if op(expr) depends on 'ω'   */
                                    /*          case of ε(t) for random var.    */
                                    if output=false then
                                            output:israndom(op(expr))
                                    )                               
                                else
                                    output:false
                            ),
                        output
                    )$
/* [wxMaxima: input   end   ] */

/* [wxMaxima: input   start ] */
ProductDecomposition(u):=
/* 
    When possible, split u into NonRandomPart / RandomPart.
    Return ** false ** if 
    ex:
         ProductDecomposition(a*K)  => false 
*/
block(
        [NonRandomPart,RandomPart,decompo],
        matchdeclare(
                NonRandomPart,
                    lambda([e], e#0 and not israndom(e)),
                RandomPart,
                    lambda([e], e#0 and israndom(e))
        ),
        defmatch(TestProdRandom,NonRandomPart*RandomPart),
        isprod:TestProdRandom(u),
        if isprod=false then
            false
        else
            [NonRandomPart,RandomPart]
    )$
/* [wxMaxima: input   end   ] */

/* [wxMaxima: input   start ] */
/*  Utilitaire pour E() */
printb(db,txt):=if db then print(txt) $
/*  Implementation of E()   */
E(u):=block(
    /*  Local variables */
    [db],
    db:false,    /* Flag pour le debug, should be VERBOSE_MODE */
    if false then u:expand(u),  /*  pour test, à éliminer.. */
    
    if israndom(u) then
        if atom(u) then
            'E(u)
        else    
            block(
            /*  Variable local  */
            [opu,argu,opdiff,opE,prodDecomp,a,foo,x],             
            /*  Initialization of local variables   */
            opu:op(u),
            argu:args(u),
            opdiff:op('diff(foo(x),x)),
            opE:op( 'E(a) ),
            /*  Try to decompose in product "" NonRandomVar * RandomVar ""    */
            prodDecomp:ProductDecomposition(u),
            if false then print("prodDecomp",prodDecomp),
            
            /************************************************************************/
            /*              Algebra rules for Expectation operator                  */            
            /************************************************************************/
            /*
             *  There are three rules that defines the Expectation operator E():
             *    (1) Idempotence EoE=E
             *    (2) Linearity for field of non-random variable
             *       a) E(X+Y)=E(X)+E(Y)
             *       b) E(aX) =a E(X)    when a is non random
             *    (3) Commutativity with partial derivative in time/space
             *
             ************************************************************************/

            if opu=opE then                
                /*  (1) Idempotence                                 */
                    E(argu[1])                        
            elseif opu="+" then
                /*  (2) Linearity                                   */
                /*      (a) Distribution with addition              */
                    map(E,u)                            
            elseif prodDecomp#false then
                /*      (b) Distribution with product of non random variable    */
                (
                printb(db,"-- Passage prodDecomp "),
                if prodDecomp[1]#1 then
                    prodDecomp[1]* E( prodDecomp[2] )
                else
                    'E( prodDecomp[2] )                        
                )
            elseif opu=opdiff  then
                /*  (3) Commutativity with derivative operator time/space       */
                (
                    printb(db,"passage d/dx or d/dt"),
                    diff( E(argu[1]) ,   argu[2],argu[3] )                            
                )
            else
            /*  Case where no simplification are encountered   */
                (printb(db,u),
                'E(u)
                )                
            )/* End of block israndom case */        
    else
        /*  No random variable  */
        u        
)$
/* [wxMaxima: input   end   ] */

/* [wxMaxima: section start ]
Validation
   [wxMaxima: section end   ] */

/* [wxMaxima: input   start ] */
if false then
    (
    /*  Fonction aléatoire */
    depends([e,ε],[x,t,ω]),
    /*  Fonction déterministe */
    depends([σ,V,ν,g,K],[x,t])
    /*declare(κ,constant)-->  no more usefull since operator E is now well defined */
    )$
/* [wxMaxima: input   end   ] */

/* [wxMaxima: input   start ] */
/**************************************  
    Validation tests for israndom       
 **************************************/
if false then validationtest([
        is(israndom(ω)=true),
        is(israndom(f(ε))=true),
        is(israndom(diff(ε,x))=true),
        is(israndom(a)=false),
        is(israndom(a(t))=false),
        is(israndom(a(t,ω))=true),   
        is(israndom(a*b)=false),     
        is(israndom(a*b(ω))=true), 
        is(israndom(ε)=true),
        /*  test E(ε) (qui n'est pas aléatoire) */
        is(israndom(E(ε))=false),
        is(israndom(diff(E(ε),x))=false),
        is(israndom(e*E(ε))=true),
        /*  Expressions compliquées */
        is(israndom(ε/(a*ε+b*diff(ε,x))*ε)=true),
        is(israndom(1/(2+1/ε))=true)
    ],"israndom")$
/* [wxMaxima: input   end   ] */

/* [wxMaxima: input   start ] */
/**************************************  
    Validation tests for ProductDecomposition       
 **************************************/
if false then validationtest([
        is(ProductDecomposition(ε/σ*1/diff(ε,x))=[1/σ,ε/diff(ε,x)]),
        is(ProductDecomposition(ε*a)=[a,ε]),
        is(ProductDecomposition(ε*E(ε))=[E(ε),ε]),
        is(ProductDecomposition(ε)=false),
        is(ProductDecomposition(a)=false),
        /*  Test with E()   */
        is(ProductDecomposition(ε*E(e))=[E(e),ε]),
        is(ProductDecomposition(ε*diff(E(ε),x)*diff(E(e**2),x))=
                            [diff(E(ε),x)*diff(E(e**2),x),ε])
    ],"ProductDecoposition")$
/* [wxMaxima: input   end   ] */

/* [wxMaxima: input   start ] */
/**************************************  
    Validation tests for E()
 **************************************/
if false then validationtest(
    [
        is(E(diff(σ,x)*32*κ) = diff(σ,x)*32*κ),
        is(map(E,[1,  e  , diff(σ,x),    e**2 , diff(e,x)**2])=
                 [1,'E(e), diff(σ,x), 'E(e**2), 'E(diff(e,x)**2)]),
        is(map(E,[  diff(e,x)   , diff(e,x,3)])=
                 [diff('E(e),x) , diff('E(e),x,3) ]),
        /*  Sommes élémentaire                  */
        is(map(E,[  σ+e  , σ+diff(e,x)-e**2       ])=
                 [σ+'E(e), σ+diff('E(e),x)-'E(e**2)]),
        /*  Multiplication élémentaire          */
        is(map(E,[σ*e**2]) = [ σ*'E(e**2)]),
        /*  Division élémentaire                */
        is(map(E,[e**2/σ,e**2*diff(e,x,5)/σ,σ/e**2,e**2/diff(e,x)])=
            ['E(e**2)/σ, 'E(e**2*diff(e,x,5))/σ, σ*'E(1/e**2), 'E(e**2/diff(e,x))]),
        /*  Validation commutation diff et /    */
        is(map(E,['diff(e**2,x)    , diff(e**2,x)   , 'diff(e**2,t)])= 
                 [diff('E(e**2),x) ,2*'E(e*diff(e,x)) , diff('E(e**2),t)]),
        /*  Validation of not full expand       */
        /*is(E((ε/diff(ε,x))*(a+ε) )=a*'E((ε/diff(ε,x)))+'E((ε**2/diff(ε,x)))),*/
        is( E((e/diff(e,x))*(a+e) )='E((e/diff(e,x))*(a+e))  ),
        /*  Validation idempotence E(E())=E()   */
        is( E( E(e**2*diff(e**2,x)) )= E(e**2*diff(e**2,x))  )        
    ],
    "E")$
/* [wxMaxima: input   end   ] */

/* Maxima can't load/batch files which end with a comment! */
"Created with wxMaxima"$
