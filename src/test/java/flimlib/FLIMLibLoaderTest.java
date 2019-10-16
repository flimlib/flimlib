/*
 * #%L
 * FLIMLib package for exponential curve fitting of fluorescence lifetime data.
 * %%
 * Copyright (C) 2010 - 2014 Gray Institute University of Oxford and Board of
 * Regents of the University of Wisconsin-Madison.
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */

package flimlib;

import java.io.IOException;

import org.junit.Test;
import org.scijava.nativelib.NativeLoader;

/**
 * Tests {@link NativeLoader}.
 * 
 * @author Dasong Gao
 */
public class FLIMLibLoaderTest {
	/** Tests {@link NativeLoader#loadLibrary}. */
	@Test
	public void testLoader() {
		try {
			NativeLoader.loadLibrary("flimlib");
			NativeLoader.loadLibrary("flimlib-jni");
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
