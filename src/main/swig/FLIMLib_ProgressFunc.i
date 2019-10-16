%module FLIMLib

%ignore t_spa_prog;
%ignore update_SPA_progress;

%inline %{
static thread_local float t_spa_prog = 0;
static void update_SPA_progress(float prog) {
	t_spa_prog = prog;
}

float getSPAProgress() {
	return t_spa_prog;
}
%}

// numinputs=0 ignores the jni argument.
%typemap(in, numinputs=0) void (*progressfunc)(float) {
	$1 = &update_SPA_progress;
}
