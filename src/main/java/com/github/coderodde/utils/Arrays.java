package com.github.coderodde.utils;

/**
 * This class exposes two radix sort implementations: one for {@code int}
 * arrays, and another one for {@code long} arrays.
 * 
 * @author Rodion "rodde" Efremov
 */
public final class Arrays {
    
    private Arrays() {}
    
    public static final LongArrays Long = new LongArrays();
    public static final IntArrays  Int  = new IntArrays();
}
