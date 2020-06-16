package com.github.coderodde.utils;

import java.util.Random;

/**
 *
 * @author Rodion "rodde" Efremov
 * @version 1.6 (Jun 13, 2020)
 */
public final class ParallelMSDRadixSortDemo {
    
    private static final int LENGTH = 50_000_000;
    static final int BUCKET_MASK = 255;
    static final int BITS_PER_BUCKET = 8;
    
    private static final int getBucketTop(final long key) {
        final int bucketIndex = (int)(key >>> 56);
        // Flip the most significant (8th) bit in order to push the negative
        // keys before the positive ones:
        return bucketIndex ^ 0b1000_0000;
    }

    private static final int getBucket(final long key, 
                                       final int recursionDepth) {
        final int bitShift = 64 - (recursionDepth + 1) * BITS_PER_BUCKET;
        return (int)(key >>> bitShift) & BUCKET_MASK;
    }
    
    public static void main(String[] args) {
//        System.out.println(Integer.toHexString(getBucketTop(0x7f00000ab0000000L)));
//        System.out.println(Integer.toHexString(getBucket(0xAABBCCL, 7)));
        runOnLongArrays(false);
        System.out.println();
        runOnLongArrays(true);
        
        long startTime = System.nanoTime();
        System.gc();
        long endTime = System.nanoTime();
        System.out.println(
                "GC in " + ((endTime - startTime) / 1000_000L) + " ms.");
        
        System.out.println("\n\n");
        
        runOnIntArrays(false);
        System.out.println();
        runOnIntArrays(true);
    }
    
    private static void runOnIntArrays(boolean printStatistics) {
        long seed = System.currentTimeMillis();
        Random random = new Random(seed);
        System.out.println("seed = " +  seed);
        
        int[] arr1 = getRandomIntArray(LENGTH, random);
        int[] arr2 = arr1.clone();
        int[] arr3 = arr1.clone();
        
        if (printStatistics) {
            System.out.println("Benchmarking on int arrays...");
        } else {
            System.out.println("Warming up on int arrays...");
        }
        
        long startTime = System.nanoTime();
        com.github.coderodde.utils.Arrays.Int.parallelSort(arr1);
        long endTime = System.nanoTime();
        
        long duration1 = (endTime - startTime) / 1_000_000L;
        
        startTime = System.nanoTime();
        java.util.Arrays.parallelSort(arr2);
        endTime = System.nanoTime();
        
        long duration2 = (endTime - startTime) / 1_000_000L;
        
        startTime = System.nanoTime();
        java.util.Arrays.sort(arr3);
        endTime = System.nanoTime();
    
        long duration3 = (endTime - startTime) / 1_000_000L;
        
        if (printStatistics) {
            System.out.println("coderodde's Arrays.parallelSort: " + duration1);
            System.out.println("Java's Arrays.parallelSort     : " + duration2);
            System.out.println("Java's Arrays.sort             : " + duration3);
            System.out.println(
                    com.github.coderodde.utils.Arrays.Int.areEqual(arr1,
                                                                   arr2,
                                                                   arr3));
        }
    }
    
    private static void runOnLongArrays(boolean printStatistics) {
        long seed = 1592224597337L; //System.currentTimeMillis();
        Random random = new Random(seed);
        System.out.println("seed = " +  seed);
        
        long[] arr1 = getRandomLongArray(LENGTH, random);
        long[] arr2 = arr1.clone();
        long[] arr3 = arr1.clone();
        
        if (printStatistics) {
            System.out.println("Benchmarking on long arrays...");
        } else {
            System.out.println("Warming up on long arrays...");
        }
        
        long startTime = System.nanoTime();
        com.github.coderodde.utils.Arrays.Long.parallelSort(arr1);
        long endTime = System.nanoTime();
        
        long duration1 = (endTime - startTime) / 1_000_000L;
        
        startTime = System.nanoTime();
        java.util.Arrays.parallelSort(arr2);
        endTime = System.nanoTime();
        
        long duration2 = (endTime - startTime) / 1_000_000L;
        
        startTime = System.nanoTime();
        java.util.Arrays.sort(arr3);
        endTime = System.nanoTime();
    
        long duration3 = (endTime - startTime) / 1_000_000L;
        
        if (printStatistics) {
            System.out.println("coderodde's Arrays.parallelSort: " + duration1);
            System.out.println("Java's Arrays.parallelSort     : " + duration2);
            System.out.println("Java's Arrays.sort             : " + duration3);
            System.out.println(
                    com.github.coderodde.utils.Arrays.Long.areEqual(arr1,
                                                                    arr2,
                                                                    arr3));
        }
    }
    
    private static long[] getRandomLongArray(int length, Random random) {
        long[] arr = new long[length];
        
        for (int i = 0; i < length; i++) {
            arr[i] = random.nextLong();
        }
        
        return arr;
    }
    
    private static int[] getRandomIntArray(int length, Random random) {
        int[] arr = new int[length];
        
        for (int i = 0; i < length; i++) {
            arr[i] = random.nextInt();
        }
        
        return arr;
    }
}
