package loci.slim;

import org.junit.Test;
import static org.junit.Assert.*;


/**
 * Tests {@link SLIMCurve} wrapper calls using swig.
 * Only tests if the calls are done correctly. Does not
 * verify that the functions work as expected. 
 * 
 * @author Zach Petersen
 */
public class SwigWrapperTest {


	final int NEVER_CALLED = -100;
	final float xInc = 0.048828125f; 
	final float y[] = {
		40.000000f, 45.000000f, 34.000000f, 54.000000f, 44.000000f,
		53.000000f, 47.000000f, 56.000000f, 62.000000f, 66.000000f,
		82.000000f, 90.000000f, 108.000000f, 122.000000f, 323.000000f,
		1155.000000f, 4072.000000f, 8278.000000f, 11919.000000f, 13152.000000f,
		13071.000000f, 12654.000000f, 11946.000000f, 11299.000000f, 10859.000000f,
		10618.000000f, 10045.000000f, 9576.000000f, 9208.000000f, 9113.000000f,
		8631.000000f, 8455.000000f, 8143.000000f, 8102.000000f, 7672.000000f,
		7384.000000f, 7463.000000f, 7254.000000f, 6980.000000f, 6910.000000f,
		6411.000000f, 6355.000000f, 6083.000000f, 5894.000000f, 5880.000000f,
		5735.000000f, 5528.000000f, 5343.000000f, 5224.000000f, 4933.000000f,
		5026.000000f, 4914.000000f, 4845.000000f, 4681.000000f, 4426.000000f,
		4485.000000f, 4271.000000f, 4295.000000f, 4183.000000f, 3989.000000f,
		3904.000000f, 3854.000000f, 3801.000000f, 3600.000000f, 3595.000000f,
		3434.000000f, 3457.000000f, 3291.000000f, 3280.000000f, 3178.000000f,
		3132.000000f, 2976.000000f, 2973.000000f, 2940.000000f, 2770.000000f,
		2969.000000f, 2851.000000f, 2702.000000f, 2677.000000f, 2460.000000f,
		2536.000000f, 2528.000000f, 2347.000000f, 2382.000000f, 2380.000000f,
		2234.000000f, 2251.000000f, 2208.000000f, 2115.000000f, 2136.000000f,
		2000.000000f, 2006.000000f, 1970.000000f, 1985.000000f, 1886.000000f,
		1898.000000f, 1884.000000f, 1744.000000f, 1751.000000f, 1797.000000f,
		1702.000000f, 1637.000000f, 1547.000000f, 1526.000000f, 1570.000000f,
		1602.000000f, 1557.000000f, 1521.000000f, 1417.000000f, 1391.000000f,
		1332.000000f, 1334.000000f, 1290.000000f, 1336.000000f, 1297.000000f,
		1176.000000f, 1189.000000f, 1220.000000f, 1209.000000f, 1217.000000f,
		1140.000000f, 1079.000000f, 1059.000000f, 1074.000000f, 1061.000000f,
		1013.000000f, 1075.000000f, 1021.000000f, 1012.000000f, 940.000000f,
		982.000000f, 866.000000f, 881.000000f, 901.000000f, 883.000000f,
		893.000000f, 845.000000f, 819.000000f, 831.000000f, 758.000000f,
		794.000000f, 779.000000f, 772.000000f, 779.000000f, 791.000000f,
		729.000000f, 732.000000f, 687.000000f, 690.000000f, 698.000000f,
		661.000000f, 647.000000f, 668.000000f, 642.000000f, 619.000000f,
		629.000000f, 656.000000f, 579.000000f, 579.000000f, 600.000000f,
		563.000000f, 584.000000f, 531.000000f, 554.000000f, 526.000000f,
		484.000000f, 530.000000f, 515.000000f, 493.000000f, 502.000000f,
		479.000000f, 445.000000f, 439.000000f, 466.000000f, 431.000000f,
		423.000000f, 451.000000f, 412.000000f, 415.000000f, 393.000000f,
		404.000000f, 390.000000f, 398.000000f, 352.000000f, 394.000000f,
		376.000000f, 338.000000f, 377.000000f, 367.000000f, 355.000000f,
		352.000000f, 375.000000f, 339.000000f, 347.000000f, 316.000000f,
		295.000000f, 322.000000f, 311.000000f, 294.000000f, 304.000000f,
		264.000000f, 293.000000f, 294.000000f, 283.000000f	
	}; 
	final int fitStart = 1;
	final int fitEnd = 203;
	final float instr[] = {
			0.00911042001098394393920898437500000f,
			0.03882249817252159118652343750000000f,
			0.13171100616455078125000000000000000f,
			0.25238901376724243164062500000000000f,
			0.27722999453544616699218750000000000f,
			0.18023300170898437500000000000000000f,
			0.07927619665861129760742187500000000f,
			0.03122770041227340698242187500000000f
	};  
	final int nInstr = 8;
	final float sig[] = {}; 
	final int paramFree[] = {1, 1, 1};
	final int nParam = 2;	
	final float chisq_trans[] = {}; 
	final int ndata = 203;
	final int ntrans = -1;
	final int ftype = 0;
	final float chisq_delta = 0;
	final int drop_bad_transients = 0;
	final int gparam[] = {};
	float[] trans = {2,3,4};
	float[] param = {0, 1000, 2}; //z, a, tau
	float[] residuals = {0};
	int[] df = {0};
	float[] chisq_global = {0};
	int noise = 5;
	float[] fitted = {0};
	int restrain = 1;
	String fitFunc = "GCI_MULTIEXP_LAMBDA";
	
	@Test
	public void testGCIGlobalWrapperCall() {		
		int result = NEVER_CALLED;
		result = SLIMCurve.GCI_marquardt_global_exps_instr(xInc, trans, ndata, ntrans, fitStart, 
				fitEnd, instr, nInstr, noise, sig, ftype, param, paramFree, nParam, restrain, 
				chisq_delta, fitted, residuals, chisq_trans, chisq_global, df, drop_bad_transients);
		System.out.println("result: " + result);
		assertTrue(result != NEVER_CALLED);
		
	}
	@Test
	public void testGCIGlobalGenericCall() {
		int result = NEVER_CALLED;
		result = SLIMCurve.GCI_marquardt_global_generic_instr(xInc, trans, ndata, ntrans, fitStart, fitEnd,
				instr, nInstr, noise, sig, param, paramFree, nParam, gparam, restrain, chisq_delta, fitFunc, fitted, 
				residuals, chisq_trans, chisq_global, df);
		System.out.println("generic result: " + result);
		assertTrue(result != NEVER_CALLED);
	}
	
	@Test
	public void testPhasorWrapperCall() {
		float z = 0;
		float u[] = {0}; 
		float v[] = {0};
		float taup[] = {0};
		float taum[] = {0};
		float tau[] = {0};
		float fitted[] = {0};
		float residuals[] = {0};
		float chiSquare[] = {0};
		
		float result = SLIMCurve.GCI_Phasor(xInc, y, fitStart, fitEnd, z, u, v, taup, taum, tau, fitted, residuals, chiSquare);
		System.out.println("phasor result: " + result + ",  phasor period: " + cLibrary.GCI_Phasor_getPeriod() + " chisq: " + chiSquare[0]);
		assertTrue(result >= 0 || result <= -5);
	}
}	