package com.github.coderodde.utils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

/**
 * This class implements an MSD (most-significant digit) radix sort for 
 * {@code long} arrays.
 * 
 * @author Rodion "rodde" Efremov
 */
public final class LongArrays {
    
    LongArrays() {}
    
    /**
     * Number of bits by which this sort processes at each recursion level.
     */
    private static final int BITS_PER_BUCKET = 8;
    
    /**
     * Number of buckets over which the current values are distributed.
     */
    private static final int BUCKETS = 1 << BITS_PER_BUCKET;
    
    /**
     * The mask for extracting the (shifted) sign bit.
     */
    private static final int SIGN_BIT_MASK = 0b1000_0000;
    
    /**
     * The mask for extracting the appropriate bucket index.
     */
    private static final int BUCKET_MASK = BUCKETS - 1;
    
    /**
     * The minimum number of keys to sort for a thread.
     */
    private static final int THREAD_THRESHOLD = 65536;
    
    /**
     * The maximum number of elements for sorting via merge sort.
     */
    private static final int MERGESORT_THRESHOLD = 4096;

    public static void parallelSort(final long[] array) {
        parallelSort(array, 0, array.length);
    }

    public static void parallelSort(final long[] array,
                                    final int fromIndex,
                                    final int toIndex) {
        final int RANGE_LENGTH = toIndex - fromIndex;

        if (RANGE_LENGTH < 2) {
            return;
        }

        final long[] buffer = array.clone();
        int threads = Math.min(RANGE_LENGTH / THREAD_THRESHOLD, 
                               Runtime.getRuntime().availableProcessors());
        
        threads = Math.max(1, threads);
        
        parallelSortImplTop(array,
                            buffer, 
                            threads, 
                            fromIndex, 
                            toIndex);
    }

    public static final boolean areEqual(final long[]... arrays) {
        for (int i = 0; i < arrays.length - 1; ++i) {
            if (arrays[i].length != arrays[i + 1].length) {
                return false;
            }
        }

        for (int i = 0; i < arrays[0].length; ++i) {
            for (int j = 0; j < arrays.length - 1; ++j) {
                if (arrays[j][i] != arrays[j + 1][i]) {
                    return false;
                }
            }
        }

        return true;
    }

    public static final boolean isSorted(final long[] array, 
                                         final int fromIndex,
                                         final int toIndex) {
        
        for (int i = fromIndex; i < toIndex - 1; ++i) {
            if (array[i] > array[i + 1]) {
                return false;
            }
        }

        return true;
    }

    public static final boolean isSorted(final long[] array) {
        return isSorted(array, 0, array.length);       
    }

    private static final void sortImpl(final long[] source,
                                       final long[] target,
                                       final int recursionDepth,
                                       final int fromIndex,
                                       final int toIndex) {
        // Try merge sort.
        if (toIndex - fromIndex <= MERGESORT_THRESHOLD) {
            mergesortAndCleanUp(source, 
                                target, 
                                recursionDepth, 
                                fromIndex,
                                toIndex);
            return;
        }

        final int[] bucketSizeMap = new int[BUCKETS];
        final int[] startIndexMap = new int[BUCKETS];
        final int[] processedMap  = new int[BUCKETS];

        // Compute the size of each bucket.
        for (int i = fromIndex; i < toIndex; ++i) {
            bucketSizeMap[getBucket(source[i], recursionDepth)]++;
        }

        // Initialize the start index map.
        startIndexMap[0] = fromIndex;

        // Compute the start index map in its entirety.
        for (int i = 1; i != BUCKETS; ++i) {
            startIndexMap[i] = startIndexMap[i - 1] +
                               bucketSizeMap[i - 1];
        }

        // Insert the entries from 'source' into their respective 'target'.
        for (int i = fromIndex; i < toIndex; ++i) {
            final long key = source[i];
            final int index = getBucket(key, recursionDepth);
            target[startIndexMap[index] + processedMap[index]++] = key;
        }

        if (recursionDepth == 7) {
            // There is nowhere to recur, return.
            return;
        }

        // Recur to sort each bucket.
        for (int i = 0; i != BUCKETS; ++i) {
            if (bucketSizeMap[i] != 0) {
                sortImpl(target,
                         source,
                         recursionDepth + 1,
                         startIndexMap[i],
                         startIndexMap[i] + bucketSizeMap[i]);
            }
        }
    }
    private static final void sortImplTop(final long[] source,
                                          final long[] target,
                                          final int fromIndex,
                                          final int toIndex) {
        // Try merge sort.
        if (toIndex - fromIndex <= MERGESORT_THRESHOLD) {
            mergesortAndCleanUp(source, 
                                target, 
                                0, 
                                fromIndex,
                                toIndex);
            return;
        }

        final int[] bucketSizeMap = new int[BUCKETS];
        final int[] startIndexMap = new int[BUCKETS];
        final int[] processedMap  = new int[BUCKETS];

        // Compute the size of each bucket.
        for (int i = fromIndex; i < toIndex; ++i) {
            bucketSizeMap[getBucketTop(source[i])]++;
        }

        // Initialize the start index map.
        startIndexMap[0] = fromIndex;

        // Compute the start index map in its entirety.
        for (int i = 1; i != BUCKETS; ++i) {
            startIndexMap[i] = startIndexMap[i - 1] +
                               bucketSizeMap[i - 1];
        }

        // Insert the entries from 'source' into their respective 'target'.
        for (int i = fromIndex; i < toIndex; ++i) {
            final long key = source[i];
            final int index = getBucketTop(key);
            target[startIndexMap[index] + processedMap[index]++] = key;
        }

        // Recur to sort each bucket.
        for (int i = 0; i != BUCKETS; ++i) {
            if (bucketSizeMap[i] != 0) {
                sortImpl(target,
                         source,
                         1,
                         startIndexMap[i],
                         startIndexMap[i] + bucketSizeMap[i]);
            }
        }
    }

    private static final boolean mergesort(final long[] source,
                                           final long[] target,
                                           final int fromIndex,
                                           final int toIndex) {
        final int RANGE_LENGTH = toIndex - fromIndex;

        long[] s = source;
        long[] t = target;

        int passes = 0;

        for (int width = 1; width < RANGE_LENGTH; width <<= 1) {
            ++passes;
            int c = 0;

            for (; c < RANGE_LENGTH / width; c += 2) {
                int left = fromIndex + c * width;
                int right = left + width;
                int i = left;

                final int leftBound = right;
                final int rightBound = Math.min(toIndex, right + width);

                while (left < leftBound && right < rightBound) {
                    t[i++] = s[right] < s[left] ?
                             s[right++] :
                             s[left++];
                }
                //! System.arraycopy here please.
                while (left < leftBound)   { t[i++] = s[left++]; }
                while (right < rightBound) { t[i++] = s[right++]; }
            }

            if (c * width < RANGE_LENGTH) {
                //! System.arraycopy here please.
                for (int i = fromIndex + c * width; i < toIndex; ++i) {
                    t[i] = s[i];
                }
            }

            final long[] tmp = s;
            s = t;
            t = tmp;
        }

        return (passes & 1) == 0;
    }

    private static final 
        void mergesortAndCleanUp(final long[] source,
                                 final long[] target,
                                 final int recursionDepth,
                                 final int fromIndex, 
                                 final int toIndex) {
        final boolean evenNumberOfMergePasses = 
                mergesort(source,
                          target, 
                          fromIndex, 
                          toIndex);

        if (evenNumberOfMergePasses) {
            // source contains the sorted range.
            if ((recursionDepth & 1) == 1) {
                // source is buffer, copy to target.
                System.arraycopy(source,
                                 fromIndex, 
                                 target,
                                 fromIndex, 
                                 toIndex - fromIndex);
            }
        } else {
            // target contains the sorted range.
            if ((recursionDepth & 1) == 0) {
                // target is buffer, copy to source.
                System.arraycopy(target, 
                                 fromIndex,
                                 source, 
                                 fromIndex, 
                                 toIndex - fromIndex);
            }
        }
    }

    private static final class BucketSizeCounter extends Thread {

        int[] localBucketSizeMap;
        private final long[] source;
        private final int recursionDepth;
        private final int fromIndex;
        private final int toIndex;

        BucketSizeCounter(final long[] source,
                          final int recursionDepth,
                          final int fromIndex,
                          final int toIndex) {
            this.source = source;
            this.recursionDepth = recursionDepth;
            this.fromIndex = fromIndex;
            this.toIndex = toIndex;
        }

        @Override
        public void run() {
            this.localBucketSizeMap = new int[BUCKETS];

            for (int i = fromIndex; i < toIndex; ++i) {
                localBucketSizeMap[getBucket(source[i], recursionDepth)]++;
            }
        }
    }
    
    private static final class TopBucketSizeCounter extends Thread {

        int[] localBucketSizeMap;
        private final long[] source;
        private final int fromIndex;
        private final int toIndex;

        TopBucketSizeCounter(final long[] source,
                             final int fromIndex,
                             final int toIndex) {
            this.source = source;
            this.fromIndex = fromIndex;
            this.toIndex = toIndex;
        }

        @Override
        public void run() {
            this.localBucketSizeMap = new int[BUCKETS];

            for (int i = fromIndex; i < toIndex; ++i) {
                localBucketSizeMap[getBucketTop(source[i])]++;
            }
        }
    }

    private static final class TopBucketInserter extends Thread {

        private final int[] startIndexMap;
        private final int[] processedMap;
        private final long[] source;
        private final long[] target;
        private final int fromIndex;
        private final int toIndex;

        TopBucketInserter(final int[] startIndexMap,
                          final int[] processedMap,
                          final long[] source,
                          final long[] target,
                          final int fromIndex,
                          final int toIndex) {
            this.startIndexMap = startIndexMap;
            this.processedMap = processedMap;
            this.source = source;
            this.target = target;
            this.fromIndex = fromIndex;
            this.toIndex = toIndex;
        }

        @Override
        public void run() {
            for (int i = fromIndex; i < toIndex; ++i) {
                final long key = source[i];
                final int index = getBucketTop(key);
                target[startIndexMap[index] + processedMap[index]++] = key;
            }
        }
    }

    private static final class BucketInserter extends Thread {

        private final int[] startIndexMap;
        private final int[] processedMap;
        private final long[] source;
        private final long[] target;
        private final int recursionDepth;
        private final int fromIndex;
        private final int toIndex;

        BucketInserter(final int[] startIndexMap,
                       final int[] processedMap,
                       final long[] source,
                       final long[] target,
                       final int recursionDepth,
                       final int fromIndex,
                       final int toIndex) {
            this.startIndexMap = startIndexMap;
            this.processedMap = processedMap;
            this.source = source;
            this.target = target;
            this.recursionDepth = recursionDepth;
            this.fromIndex = fromIndex;
            this.toIndex = toIndex;
        }

        @Override
        public void run() {
            for (int i = fromIndex; i < toIndex; ++i) {
                final long key = source[i];
                final int index = getBucket(key, recursionDepth);
                target[startIndexMap[index] + processedMap[index]++] = key;
            }
        }
    }

    private static final class Sorter extends Thread {

        private final List<Task> taskList;

        Sorter(final List<Task> taskList) {
            this.taskList = taskList;
        }

        @Override
        public void run() {
            for (final Task task : taskList) {
                // Choose parallel or sequential.
                if (task.threads > 1) {
                    if (task.recursionDepth == 0) {
                        parallelSortImplTop(task.source,
                                            task.target,
                                            task.threads,
                                            task.fromIndex,
                                            task.toIndex);
                    } else {
                        parallelSortImpl(task.source,
                                         task.target,
                                         task.threads,
                                         task.recursionDepth,
                                         task.fromIndex,
                                         task.toIndex);
                    }
                } else {
                    if (task.recursionDepth == 0) {
                        sortImplTop(task.source, 
                                    task.target, 
                                    task.fromIndex, 
                                    task.toIndex);
                    } else {
                        sortImpl(task.source,
                                 task.target,
                                 task.recursionDepth,
                                 task.fromIndex,
                                 task.toIndex);
                    }
                }
            }
        }
    }

    private static final class Task {

        private final long[] source;
        private final long[] target;
        private final int threads;
        private final int recursionDepth;
        private final int fromIndex;
        private final int toIndex;

        Task(final long[] source,
             final long[] target,
             final int threads,
             final int recursionDepth,
             final int fromIndex,
             final int toIndex) {
            this.source = source;
            this.target = target;
            this.threads = threads;
            this.recursionDepth = recursionDepth;
            this.fromIndex = fromIndex;
            this.toIndex = toIndex;
        }
    }

    private static final void parallelSortImpl(final long[] source,
                                               final long[] target,
                                               final int threads,
                                               final int recursionDepth,
                                               final int fromIndex,
                                               final int toIndex) {
        final int RANGE_LENGTH = toIndex - fromIndex;

        if (RANGE_LENGTH <= MERGESORT_THRESHOLD) {
            mergesortAndCleanUp(source, 
                                target, 
                                recursionDepth, 
                                fromIndex, 
                                toIndex);
            return;
        }

        if (threads < 2) {
            sortImpl(source, target, recursionDepth, fromIndex, toIndex);
            return;
        }

        // Create the bucket size counter threads.
        final BucketSizeCounter[] counters = new BucketSizeCounter[threads];
        final int SUB_RANGE_LENGTH = RANGE_LENGTH / threads;
        int start = fromIndex;

        for (int i = 0; i != threads - 1; ++i, start += SUB_RANGE_LENGTH) {
            counters[i] = new BucketSizeCounter(source,
                                                recursionDepth,
                                                start,
                                                start + SUB_RANGE_LENGTH);
            counters[i].start();
        }

        counters[threads - 1] = 
                new BucketSizeCounter(source,
                                      recursionDepth,
                                      start,
                                      toIndex);

        // Run the last counter in this thread while other are already on their
        // way.
        counters[threads - 1].run();

        try {
            for (int i = 0; i != threads - 1; ++i) {
                counters[i].join();
            }
        } catch (final InterruptedException ie) {
            ie.printStackTrace();
            return;
        }

        final int[] bucketSizeMap = new int[BUCKETS];
        final int[] startIndexMap = new int[BUCKETS];

        // Count the size of each processed bucket.
        for (int i = 0; i != threads; ++i) {
            for (int j = 0; j != BUCKETS; ++j) {
                bucketSizeMap[j] += counters[i].localBucketSizeMap[j];
            }
        }

        // Prepare the starting indices of each bucket.
        startIndexMap[0] = fromIndex;

        for (int i = 1; i != BUCKETS; ++i) {
            startIndexMap[i] = startIndexMap[i - 1] +
                               bucketSizeMap[i - 1];
        }

        // Create the inserter threads.
        final BucketInserter[] inserters = new BucketInserter[threads - 1];
        final int[][] processedMaps = new int[threads][BUCKETS];

        // Make processedMaps of each thread independent of the other.
        for (int i = 1; i != threads; ++i) {
            int[] partialBucketSizeMap = counters[i - 1].localBucketSizeMap;

            for (int j = 0; j != BUCKETS; ++j) {
                processedMaps[i][j] = 
                        processedMaps[i - 1][j] + partialBucketSizeMap[j];
            }
        }

        int startIndex = fromIndex;

        for (int i = 0; i != threads - 1; ++i, startIndex += SUB_RANGE_LENGTH) {
            inserters[i] =
                    new BucketInserter(startIndexMap,
                                       processedMaps[i],
                                       source,
                                       target,
                                       recursionDepth,
                                       startIndex,
                                       startIndex + SUB_RANGE_LENGTH);
            inserters[i].start();
        }

        // Run the last inserter in this thread while other are on their ways.
        new BucketInserter(startIndexMap,
                           processedMaps[threads - 1],
                           source,
                           target,
                           recursionDepth,
                           startIndex,
                           toIndex).run();

        try {
            for (int i = 0; i != threads - 1; ++i) {
                inserters[i].join();
            }
        } catch (final InterruptedException ie) {
            ie.printStackTrace();
            return;
        }

        if (recursionDepth == 7) {
            // Nowhere to recur.
            return;
        }

        int nonEmptyBucketAmount = 0;

        for (int i : bucketSizeMap) {
            if (i != 0) {
                ++nonEmptyBucketAmount;
            }
        }

        final int SPAWN_DEGREE = Math.min(nonEmptyBucketAmount, threads);
        final List<Integer>[] bucketIndexListArray = new List[SPAWN_DEGREE];

        for (int i = 0; i != SPAWN_DEGREE; ++i) {
            bucketIndexListArray[i] = new ArrayList<>(nonEmptyBucketAmount);
        }

        final int[] threadCountMap = new int[SPAWN_DEGREE];

        for (int i = 0; i != SPAWN_DEGREE; ++i) {
            threadCountMap[i] = threads / SPAWN_DEGREE;
        }

        for (int i = 0; i != threads % SPAWN_DEGREE; ++i) {
            ++threadCountMap[i];
        }

        final List<Integer> nonEmptyBucketIndices = 
                new ArrayList<>(nonEmptyBucketAmount);


        for (int i = 0; i != BUCKETS; ++i) {
            if (bucketSizeMap[i] != 0) {
                nonEmptyBucketIndices.add(i);
            }
        }

        Collections.sort(nonEmptyBucketIndices, 
                         new BucketSizeComparator(bucketSizeMap));

        final int OPTIMAL_SUBRANGE_LENGTH = RANGE_LENGTH / SPAWN_DEGREE;
        int listIndex = 0;
        int packed = 0;
        int f = 0;
        int j = 0;

        while (j < nonEmptyBucketIndices.size()) {
            int tmp = bucketSizeMap[nonEmptyBucketIndices.get(j++)];
            packed += tmp;

            if (packed >= OPTIMAL_SUBRANGE_LENGTH
                    || j == nonEmptyBucketIndices.size()) {
                packed = 0;

                for (int i = f; i < j; ++i) {
                    bucketIndexListArray[listIndex]
                            .add(nonEmptyBucketIndices.get(i));
                }

                ++listIndex;
                f = j;
            }
        }

        final Sorter[] sorters = new Sorter[SPAWN_DEGREE];
        final List<List<Task>> llt = new ArrayList<>(SPAWN_DEGREE);

        for (int i = 0; i != SPAWN_DEGREE; ++i) {
            final List<Task> lt = new ArrayList<>();

            for (int idx : bucketIndexListArray[i]) {
                lt.add(new Task(target,
                                source,
                                threadCountMap[i],
                                recursionDepth + 1,
                                startIndexMap[idx],
                                startIndexMap[idx] + bucketSizeMap[idx]));
            }

            llt.add(lt);
        }

        for (int i = 0; i != SPAWN_DEGREE - 1; ++i) {
            sorters[i] = new Sorter(llt.get(i));
            sorters[i].start();
        }

        new Sorter(llt.get(SPAWN_DEGREE - 1)).run();

        try {
            for (int i = 0; i != SPAWN_DEGREE - 1; ++i) {
                sorters[i].join();
            }
        } catch (final InterruptedException ie) {
            ie.printStackTrace();
            return;
        }
    }

    private static final void parallelSortImplTop(final long[] source,
                                                  final long[] target,
                                                  final int threads,
                                                  final int fromIndex,
                                                  final int toIndex) {
        final int RANGE_LENGTH = toIndex - fromIndex;

        if (RANGE_LENGTH <= MERGESORT_THRESHOLD) {
            mergesortAndCleanUp(source, 
                                target, 
                                0, 
                                fromIndex, 
                                toIndex);
            return;
        }

        if (threads < 2) {
            sortImplTop(source, 
                        target, 
                        fromIndex, 
                        toIndex);
            return;
        }

        // Create the bucket size counter threads.
        final TopBucketSizeCounter[] counters = 
                new TopBucketSizeCounter[threads];
        
        final int SUB_RANGE_LENGTH = RANGE_LENGTH / threads;
        int start = fromIndex;

        for (int i = 0; i != threads - 1; ++i, start += SUB_RANGE_LENGTH) {
            counters[i] = new TopBucketSizeCounter(source,
                                                   start,
                                                   start + SUB_RANGE_LENGTH);
            counters[i].start();
        }

        counters[threads - 1] = 
                new TopBucketSizeCounter(source,
                                         start,
                                         toIndex);

        // Run the last counter in this thread while other are already on their
        // way.
        counters[threads - 1].run();

        try {
            for (int i = 0; i != threads - 1; ++i) {
                counters[i].join();
            }
        } catch (final InterruptedException ie) {
            ie.printStackTrace();
            return;
        }

        final int[] bucketSizeMap = new int[BUCKETS];
        final int[] startIndexMap = new int[BUCKETS];

        // Count the size of each processed bucket.
        for (int i = 0; i != threads; ++i) {
            for (int j = 0; j != BUCKETS; ++j) {
                bucketSizeMap[j] += counters[i].localBucketSizeMap[j];
            }
        }

        // Prepare the starting indices of each bucket.
        startIndexMap[0] = fromIndex;

        for (int i = 1; i != BUCKETS; ++i) {
            startIndexMap[i] = startIndexMap[i - 1] +
                               bucketSizeMap[i - 1];
        }

        // Create the inserter threads.
        final TopBucketInserter[] inserters = 
                new TopBucketInserter[threads - 1];
        
        final int[][] processedMaps = new int[threads][BUCKETS];

        // Make processedMaps of each thread independent of the other.
        for (int i = 1; i != threads; ++i) {
            int[] partialBucketSizeMap = counters[i - 1].localBucketSizeMap;

            for (int j = 0; j != BUCKETS; ++j) {
                processedMaps[i][j] = 
                        processedMaps[i - 1][j] + partialBucketSizeMap[j];
            }
        }

        int startIndex = fromIndex;

        for (int i = 0; i != threads - 1; ++i, startIndex += SUB_RANGE_LENGTH) {
            inserters[i] =
                    new TopBucketInserter(startIndexMap,
                                          processedMaps[i],
                                          source,
                                          target,
                                          startIndex,
                                          startIndex + SUB_RANGE_LENGTH);
            inserters[i].start();
        }

        // Run the last inserter in this thread while other are on their ways.
        new TopBucketInserter(startIndexMap,
                              processedMaps[threads - 1],
                              source,
                              target,
                              startIndex,
                              toIndex).run();

        try {
            for (int i = 0; i != threads - 1; ++i) {
                inserters[i].join();
            }
        } catch (final InterruptedException ie) {
            ie.printStackTrace();
            return;
        }

        int nonEmptyBucketAmount = 0;

        for (int i : bucketSizeMap) {
            if (i != 0) {
                ++nonEmptyBucketAmount;
            }
        }

        final int SPAWN_DEGREE = Math.min(nonEmptyBucketAmount, threads);
        final List<Integer>[] bucketIndexListArray = new List[SPAWN_DEGREE];

        for (int i = 0; i != SPAWN_DEGREE; ++i) {
            bucketIndexListArray[i] = new ArrayList<>(nonEmptyBucketAmount);
        }

        final int[] threadCountMap = new int[SPAWN_DEGREE];

        for (int i = 0; i != SPAWN_DEGREE; ++i) {
            threadCountMap[i] = threads / SPAWN_DEGREE;
        }

        for (int i = 0; i != threads % SPAWN_DEGREE; ++i) {
            ++threadCountMap[i];
        }

        final List<Integer> nonEmptyBucketIndices = 
                new ArrayList<>(nonEmptyBucketAmount);


        for (int i = 0; i != BUCKETS; ++i) {
            if (bucketSizeMap[i] != 0) {
                nonEmptyBucketIndices.add(i);
            }
        }

        Collections.sort(nonEmptyBucketIndices, 
                         new BucketSizeComparator(bucketSizeMap));

        final int OPTIMAL_SUBRANGE_LENGTH = RANGE_LENGTH / SPAWN_DEGREE;
        int listIndex = 0;
        int packed = 0;
        int f = 0;
        int j = 0;

        while (j < nonEmptyBucketIndices.size()) {
            int tmp = bucketSizeMap[nonEmptyBucketIndices.get(j++)];
            packed += tmp;

            if (packed >= OPTIMAL_SUBRANGE_LENGTH
                    || j == nonEmptyBucketIndices.size()) {
                packed = 0;

                for (int i = f; i < j; ++i) {
                    bucketIndexListArray[listIndex]
                            .add(nonEmptyBucketIndices.get(i));
                }

                ++listIndex;
                f = j;
            }
        }

        final Sorter[] sorters = new Sorter[SPAWN_DEGREE];
        final List<List<Task>> llt = new ArrayList<>(SPAWN_DEGREE);

        for (int i = 0; i != SPAWN_DEGREE; ++i) {
            final List<Task> lt = new ArrayList<>();

            for (int idx : bucketIndexListArray[i]) {
                lt.add(new Task(target,
                                source,
                                threadCountMap[i],
                                1,
                                startIndexMap[idx],
                                startIndexMap[idx] + bucketSizeMap[idx]));
            }

            llt.add(lt);
        }

        for (int i = 0; i != SPAWN_DEGREE - 1; ++i) {
            sorters[i] = new Sorter(llt.get(i));
            sorters[i].start();
        }

        new Sorter(llt.get(SPAWN_DEGREE - 1)).run();

        try {
            for (int i = 0; i != SPAWN_DEGREE - 1; ++i) {
                sorters[i].join();
            }
        } catch (final InterruptedException ie) {
            ie.printStackTrace();
            return;
        }
    }

    private static final class BucketSizeComparator 
    implements Comparator<Integer> {
        private final int[] bucketSizeMap;

        BucketSizeComparator(final int[] bucketSizeMap) {
            this.bucketSizeMap = bucketSizeMap;
        }

        @Override
        public int compare(final Integer i1, final Integer i2) {
            final int sz1 = bucketSizeMap[i1];
            final int sz2 = bucketSizeMap[i2];
            return sz2 - sz1;
        }
    }
    
    private static final int getBucketTop(final long key) {
        final int bucketIndex = (int)(key >>> 56);
        // Flip the most significant (8th) bit in order to push the negative
        // keys before the positive ones:
        return bucketIndex ^ SIGN_BIT_MASK;
    }

    private static final int getBucket(final long key, 
                                       final int recursionDepth) {
        final int bitShift = 64 - (recursionDepth + 1) * BITS_PER_BUCKET;
        return (int)(key >>> bitShift) & BUCKET_MASK;
    }
}