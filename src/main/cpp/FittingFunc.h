typedef void (*fitfunc)(float, float [], float *, float [], int);
class FittingFunction {
public:
    fitfunc func_ptr;
    virtual void func(float, float [], float *, float [], int) = 0; 
};