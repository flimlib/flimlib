package slim;

public interface FitFunc {
	// builtin fitting functions
	public static final FitFunc GCI_MULTIEXP_LAMBDA = SLIMCurveConstants.GCI_MULTIEXP_LAMBDA;
	public static final FitFunc GCI_MULTIEXP_TAU = SLIMCurveConstants.GCI_MULTIEXP_TAU;
	public static final FitFunc GCI_STRETCHEDEXP = SLIMCurveConstants.GCI_STRETCHEDEXP;

	public float fit(final float x, final float[] param, final float[] dy_dparam);
}
