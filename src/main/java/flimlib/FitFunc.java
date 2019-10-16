package flimlib;

public interface FitFunc {
	// builtin fitting functions
	public static final FitFunc GCI_MULTIEXP_LAMBDA = FLIMLibConstants.GCI_MULTIEXP_LAMBDA;
	public static final FitFunc GCI_MULTIEXP_TAU = FLIMLibConstants.GCI_MULTIEXP_TAU;
	public static final FitFunc GCI_STRETCHEDEXP = FLIMLibConstants.GCI_STRETCHEDEXP;

	public float fit(final float x, final float[] param, final float[] dy_dparam);
}
