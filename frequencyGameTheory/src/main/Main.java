package main;
import java.io.*; //BufferedWriter, FileWriter, and File
import java.util.*; //ArrayList, HashMap, and Map

/**
 *  Goals: 
 * 1. Make a simpler way to add in the chord progression, for possible experimentation with other chord progressions.
 * 2. Fix the variance score
 *     - Perform data analysis on the variance score in actual music to find optimal variance score distribution
 *     - Fix the skipping beat 1 problem
 *     - See if there is a potential problem with the startswith function when checking A vs A#, like maybe A# is being counted as A
 *     - Potentially implement a spiral method of calculating note distance?? Or maybe we just need variance without consideration of distance.
 *     - Or maybe we need a seperate note and octave variance score, and we can play with the weight of each one.
 *     - Fix predictive harmony player and whatever is going on with that a 
 * 3. Improve Harmony Score
 *     - use the fraction approximation technique
 *     - Perform data analysis on the harmony score in actual music to find optimal harmony score distribution
 */

public class Main {
    protected static Map<String, double[]> notesFreqMap;
    protected static int[] chordProgressionFreq;

    public static void main(String[] args) throws Exception {
        // pairHarmonyScore(2630, 3048); // 5/6
        // pairHarmonyScore(2630, 3740); // 32/45
        // pairHarmonyScore(180, 256); // 32/45
        // pairHarmonyScore(45, 64); // 45/64

        buildNotesFrequenciesMap();
        Player p1 = new PredictiveHarmonyPlayer();
        Player p2 = new StepwisePlayer();
        String fileName = getFileName(p1, p2);
        BufferedWriter bw = new BufferedWriter(new FileWriter(new File(fileName))); //will replace file if simulation has already been run
        ArrayList<int[]> chordProgression = new ArrayList<int[]>(); //every int[] is one chord (usually 4 notes)
        int[] bb7 = {116, 147, 175, 208}; // Bb7 -(equal temperament approx.frequencies)-> 116.54(Bb2),146.83(D3), 174.61(F3), Ab(207.65)
        int[] eb7 = {156, 196, 233, 277}; // Eb7 -(equal temperament approx.frequencies)-> Eb3 	155.56, G3	196.00, Bb3 	233.08, Db4 	277.18
        int[] cm7 = {131, 156, 196, 233}; // Cm7 -(equal temperament approx.frequencies)-> C3	130.81, Eb3 	155.56, G3	196.00, Bb3 	233.08
        int[] f7 = {175, 220, 262, 311}; // F7 -(equal temperament approx.frequencies)-> F3	174.61, A3	220.00, C4	261.63, Eb4 	311.13
        /**
         * Add the chords in the right order
         * 12 bar Bb blues progression: Bb7, Eb7, Bb7, Bb7, Eb7, Eb7, Bb7, Bb7, Cm7, F7, Bb7, F7
         */
        for (int i = 0; i < 8; i++) { //goes through chord 8 times
            chordProgression.add(bb7);
            chordProgression.add(eb7);
            chordProgression.add(bb7);
            chordProgression.add(bb7);
            chordProgression.add(eb7);
            chordProgression.add(eb7);
            chordProgression.add(bb7);
            chordProgression.add(bb7);
            chordProgression.add(cm7);
            chordProgression.add(f7);
            chordProgression.add(bb7);
            chordProgression.add(f7);
        }
        ArrayList<Integer> allPastNotes = new ArrayList<Integer>();
        ArrayList<Integer> p1Freq = new ArrayList<Integer>(); //saves player 1 notes to list
        ArrayList<Integer> p2Freq = new ArrayList<Integer>(); //saves player 2 notes to list
        double payoffSum = 0.0;
        int measureNum = -1; // counts the measures
        for (int beatNum = 0; beatNum < 96*8; beatNum++) { // subdividing by eight notes there will be 96 beats in a 12 bar blues, and then 8 repetitions of the 12 bar blues
            if (beatNum % 8 == 0) measureNum++; // increments measureNum at the start of every 8 beats --> one measure
            chordProgressionFreq = chordProgression.get(measureNum);
            int freqOne = p1.genNote(); //gets note based on p1's strategy
            int freqTwo = p2.genNote();
            p1Freq.add(freqOne);
            p2Freq.add(freqTwo);
            /**
             * weighted average for the current note: chord progression frequency is weighted at 60%
             * Each individual player's played note's frequency is weighted at 20%
             **/
            allPastNotes.add(freqOne);
            allPastNotes.add(freqTwo);
            double varianceScore = 0;
            if (beatNum != 0) {
                //calculate payoff in here
                varianceScore = calcVarianceScore(allPastNotes);
            }
            double harmonyScore = betterCalcHarmonyScore(chordProgressionFreq, freqOne, freqTwo);
            final double normalizationFactor = 1658.833002; //calculation shown in paper; uses raw variance/harmony scores for normalization
            varianceScore*=normalizationFactor;
            harmonyScore/=3; // If you divide by 3 it sounds fire lowkey
            if (beatNum == 0) {
                varianceScore = harmonyScore; //makes payoff 0 instead of negative for the first beat
            }
            double payoff = (varianceScore - harmonyScore)/(varianceScore + harmonyScore);
            if (p1 instanceof SimpleReinforcementPlayer || p1 instanceof StepwisePlayer || p1 instanceof ChordFollowingReinforcementLearning || p1 instanceof ChordSpecificReinforcementPlayer) {
                p1.update(payoff);
            }
            if (p2 instanceof SimpleReinforcementPlayer || p2 instanceof StepwisePlayer || p2 instanceof ChordFollowingReinforcementLearning || p2 instanceof ChordSpecificReinforcementPlayer) {
                p2.update(payoff);
            }
            if (p1 instanceof DoubleReinforcementPlayer) {
                p1.update(payoff, freqTwo);
            }
            if (p2 instanceof DoubleReinforcementPlayer) {
                p2.update(payoff, freqOne);
            }
            if (p1 instanceof PredictiveHarmonyPlayer) {
                p1.update(freqTwo);
            }
            if (p2 instanceof PredictiveHarmonyPlayer) {
                p2.update(freqOne);
            }
            bw.write(varianceScore + "\t" + harmonyScore + "\t" + payoff + "\n");
            payoffSum+=payoff;
            bw.flush();
        }
        Musician.play(p1Freq, p2Freq, fileName); //generates MIDI files for player 1
        // Musician.play(p2Freq, 2, fileName); //generates MIDI files for player 2
        bw.write("\nAverage Payoff: " + String.format("%.4f", payoffSum/(96.0*8)));
        bw.flush(); //write to relevant notepad
        bw.close(); //prevent resource leaks
    }

    private static void buildNotesFrequenciesMap() {
        notesFreqMap = new LinkedHashMap<String, double[]>();
        notesFreqMap.put("A0", new double[]{28, 28.31}); //starts at 28 Hz
        notesFreqMap.put("A#0/Bb0", new double[]{28.31, 29.99});
        notesFreqMap.put("B0", new double[]{29.99, 31.77});
        notesFreqMap.put("C1", new double[]{31.77, 33.66});
        notesFreqMap.put("C#1/Db1", new double[]{33.66, 35.67});
        notesFreqMap.put("D1", new double[]{35.67, 37.78});
        notesFreqMap.put("D#1/Eb1", new double[]{37.78, 40.03});
        notesFreqMap.put("E1", new double[]{40.03, 42.41});
        notesFreqMap.put("F1", new double[]{42.41, 44.93});
        notesFreqMap.put("F#1/Gb1", new double[]{44.93, 47.61});
        notesFreqMap.put("G1", new double[]{47.61, 50.43});
        notesFreqMap.put("G#1/Ab1", new double[]{50.43, 53.43});
        notesFreqMap.put("A1", new double[]{53.43, 56.61});
        notesFreqMap.put("A#1/Bb1", new double[]{56.61, 60.00});
        notesFreqMap.put("B1", new double[]{60.00, 63.55});
        notesFreqMap.put("C2", new double[]{63.55, 67.32});
        notesFreqMap.put("C#2/Db2", new double[]{67.32, 71.34});
        notesFreqMap.put("D2", new double[]{71.34, 75.57});
        notesFreqMap.put("D#2/Eb2", new double[]{75.57, 80.07});
        notesFreqMap.put("E2", new double[]{80.07, 84.83});
        notesFreqMap.put("F2", new double[]{84.83, 89.87});
        notesFreqMap.put("F#2/Gb2", new double[]{89.87, 95.22});
        notesFreqMap.put("G2", new double[]{95.22, 100.87});
        notesFreqMap.put("G#2/Ab2", new double[]{100.87, 106.87});
        notesFreqMap.put("A2", new double[]{106.87, 113.23});
        notesFreqMap.put("A#2/Bb2", new double[]{113.23, 120.00});
        notesFreqMap.put("B2", new double[]{120.00, 127.10});
        notesFreqMap.put("C3", new double[]{127.10, 134.64});
        notesFreqMap.put("C#3/Db3", new double[]{134.64, 142.68});
        notesFreqMap.put("D3", new double[]{142.68, 151.14});
        notesFreqMap.put("D#3/Eb3", new double[]{151.14, 160.14});
        notesFreqMap.put("E3", new double[]{160.14, 169.65});
        notesFreqMap.put("F3", new double[]{169.65, 179.74});
        notesFreqMap.put("F#3/Gb3", new double[]{179.74, 190.44});
        notesFreqMap.put("G3", new double[]{190.44, 201.74});
        notesFreqMap.put("G#3/Ab3", new double[]{201.74, 213.74});
        notesFreqMap.put("A3", new double[]{213.74, 226.46});
        notesFreqMap.put("A#3/Bb3", new double[]{226.46, 240.00});
        notesFreqMap.put("B3", new double[]{240.00, 254.20});
        notesFreqMap.put("C4", new double[]{254.20, 269.29});
        notesFreqMap.put("C#4/Db4", new double[]{269.29, 285.36});
        notesFreqMap.put("D4", new double[]{285.36, 302.28});
        notesFreqMap.put("D#4/Eb4", new double[]{302.28, 320.29});
        notesFreqMap.put("E4", new double[]{320.29, 339.30});
        notesFreqMap.put("F4", new double[]{339.30, 359.48});
        notesFreqMap.put("F#4/Gb4", new double[]{359.48, 380.89});
        notesFreqMap.put("G4", new double[]{380.89, 403.49});
        notesFreqMap.put("G#4/Ab4", new double[]{403.49, 427.48});
        notesFreqMap.put("A4", new double[]{427.48, 452.92});
        notesFreqMap.put("A#4/Bb4", new double[]{452.92, 480.00});
        notesFreqMap.put("B4", new double[]{480.00, 508.40});
        notesFreqMap.put("C5", new double[]{508.40, 538.58});
        notesFreqMap.put("C#5/Db5", new double[]{538.58, 570.72});
        notesFreqMap.put("D5", new double[]{570.72, 604.56});
        notesFreqMap.put("D#5/Eb5", new double[]{604.56, 640.57});
        notesFreqMap.put("E5", new double[]{640.57, 678.60});
        notesFreqMap.put("F5", new double[]{678.60, 718.96});
        notesFreqMap.put("F#5/Gb5", new double[]{718.96, 761.77});
        notesFreqMap.put("G5", new double[]{761.77, 806.97});
        notesFreqMap.put("G#5/Ab5", new double[]{806.97, 854.96});
        notesFreqMap.put("A5", new double[]{854.96, 905.84});
        notesFreqMap.put("A#5/Bb5", new double[]{905.84, 960.00});
        notesFreqMap.put("B5", new double[]{960.00, 1016.80});
        notesFreqMap.put("C6", new double[]{1016.80, 1077.16});
        notesFreqMap.put("C#6/Db6", new double[]{1077.16, 1141.44});
        notesFreqMap.put("D6", new double[]{1141.44, 1209.12});
        notesFreqMap.put("D#6/Eb6", new double[]{1209.12, 1281.14});
        notesFreqMap.put("E6", new double[]{1281.14, 1357.20});
        notesFreqMap.put("F6", new double[]{1357.20, 1437.92});
        notesFreqMap.put("F#6/Gb6", new double[]{1437.92, 1523.55});
        notesFreqMap.put("G6", new double[]{1523.55, 1613.94});
        notesFreqMap.put("G#6/Ab6", new double[]{1613.94, 1709.92});
        notesFreqMap.put("A6", new double[]{1709.92, 1811.68});
        notesFreqMap.put("A#6/Bb6", new double[]{1811.68, 1920.00});
        notesFreqMap.put("B6", new double[]{1920.00, 2033.60});
        notesFreqMap.put("C7", new double[]{2033.60, 2154.32});
        notesFreqMap.put("C#7/Db7", new double[]{2154.32, 2282.88});
        notesFreqMap.put("D7", new double[]{2282.88, 2418.24});
        notesFreqMap.put("D#7/Eb7", new double[]{2418.24, 2562.28});
        notesFreqMap.put("E7", new double[]{2562.28, 2714.40});
        notesFreqMap.put("F7", new double[]{2714.40, 2875.84});
        notesFreqMap.put("F#7/Gb7", new double[]{2875.84, 3047.09});
        notesFreqMap.put("G7", new double[]{3047.09, 3227.88});
        notesFreqMap.put("G#7/Ab7", new double[]{3227.88, 3419.84});
        notesFreqMap.put("A7", new double[]{3419.84, 3623.36});
        notesFreqMap.put("A#7/Bb7", new double[]{3623.36, 3840.00});
        notesFreqMap.put("B7", new double[]{3840.00, 4067.20});
        notesFreqMap.put("C8", new double[]{4067.20, 4186}); //no notes past 4186
    }

    //variance when not counting octaves as the same note:

    // private static double calcVarianceScore(ArrayList<Integer> allPastNotes) {
    //     /**
    //      * Separate past notes played into buckets of actual notes
    //      * Find frequencies
    //      * Find variance of those frequencies
    //      */
    //     //Separating past notes played into buckets of actual notes
    //     int[] noteCounts = new int[notesFreqMap.size()];
    //     for (int x: allPastNotes) {
    //         int i = 0;
    //         for (String s: notesFreqMap.keySet()) {
    //             if (x >= notesFreqMap.get(s)[0] && x <= notesFreqMap.get(s)[1]) {
    //                 noteCounts[i]++;
    //                 break;
    //             }
    //             i++;
    //         }
    //     }
    //     //Calculating Shannon Diversity Index
    //     double diversityIndex = 0;
    //     for (int x: noteCounts) {
    //         double Pi = (double)x / allPastNotes.size();
    //         double calc = -Pi*Math.log(Pi); //-Pi * ln(Pi)
    //         if (Double.isNaN(calc)) {
    //             calc = 0;
    //         }
    //         diversityIndex+=calc;
    //     }
    //     return diversityIndex/Math.log(noteCounts.length); //H/Hmax
    // }

    //counts octaves as the same note
    private static double calcVarianceScore(ArrayList<Integer> allPastNotes) {
        Map<String, Integer> noteCounts = new HashMap<String, Integer>();
        final String[] notes = new String[]{"A", "A#", "B", "C", "C#", "D", "D#", "E", "F", "F#", "G", "G#"};
        for (String s: notes) {
            noteCounts.put(s, 0);
        }
        for (int x: allPastNotes) {
            for (String s: notesFreqMap.keySet()) {
                if (x >= notesFreqMap.get(s)[0] && x < notesFreqMap.get(s)[1]) {
                    if (s.contains("#")) {
                        for (int i =notes.length-1; i >=0; i--) {
                            String note = notes[i];
                            if (s.startsWith(note)) {
                                noteCounts.put(note, noteCounts.get(note)+1);
                                // System.out.println("Actuall note: " + s);
                                // System.out.println("Recorded Note: " + note);
                                break;
                            }
                        }
                    }
                    else {
                        for (String note : notes) {
                            if (s.startsWith(note)) {
                                noteCounts.put(note, noteCounts.get(note)+1);
                                // System.out.println("Actual note: " + s);
                                // System.out.println("Recorded Note: " + note);
                                break;
                            }
                        }
                    }
                }
            }
        }
        //Calculating Shannon Diversity Index
        double diversityIndex = 0;
        for (int x: noteCounts.values()) {
            double Pi = (double)x / allPastNotes.size();
            double calc = -Pi*Math.log(Pi); //-Pi * ln(Pi)
            if (Double.isNaN(calc)) {
                calc = 0;
            }
            diversityIndex+=calc;
        }
        return diversityIndex/Math.log(noteCounts.size()); //H/Hmax
    }

    private static double calcHarmonyScore(int[] chord, int freqOne, int freqTwo) {
        ArrayList<Integer> allNotes = new ArrayList<Integer>();
        for (int x: chord) {
            allNotes.add(x);
        }
        allNotes.add(freqOne); allNotes.add(freqTwo);
        int sum = 0;
        for (int i = 0; i < allNotes.size()-1; i++) {
            for (int j = i + 1; j < allNotes.size(); j++) {
                sum+=addNumDenomSimplifiedFraction(allNotes.get(i), allNotes.get(j));
            }
        }
        int numTimesOfLoop = (allNotes.size()*(allNotes.size()-1))/2;
        return (int) (Math.round(sum/((double)(numTimesOfLoop))));
    }

    /**
     * Takes in two integers, x and y
     * Simplifies the fraction x/y and returns the sum of the numerator and denominator
     */
    private static int addNumDenomSimplifiedFraction(int x, int y) {
        int gcd = findGCD(x,y);
        return (x/gcd) + (y/gcd);
    }

    //Credit for Algorithm: Geeks for Geeks
    //https://www.geeksforgeeks.org/program-to-find-gcd-or-hcf-of-two-numbers/
    private static int findGCD(int x, int y) {
        int result = Math.min(x, y);
        while (result > 0) {
            if (x % result == 0 && y % result == 0) {
                break;
            }
            result--;
        }
        return result;
    }

    private static int betterCalcHarmonyScore(int[] chord, int freqOne, int freqTwo) {
        int[] allNotes = new int[chord.length+2];
        for (int i = 0; i < chord.length; i++) {
            allNotes[i] = chord[i];
        }
        allNotes[chord.length] = freqOne;
        allNotes[chord.length+1] = freqTwo;
        int sum = 0;
        for (int i = 0; i < allNotes.length-1; i++) {
            for (int j = i + 1; j < allNotes.length; j++) {
                sum+=pairHarmonyScore(allNotes[i], allNotes[j]);
            }
        }
        return sum;
    }

    /**
     * Takes in two frequencies, freqOne and freqTwo, and finds the simplest fraction
     * of the two frequencies with a 1% error margin.
     * @param freqOne
     * @param freqTwo
     * @return Sum of the numerator and the demoninator.
     */
    private static int pairHarmonyScore(int freqOne, int freqTwo) {
        int numerator; int denominator;
        // Let's say that frequency1 = 200, and freq2 = 401. So we have the fraction 200/401. 
        // Freq1 must always be less than freq2
        if (freqOne > freqTwo) {
            int temp = freqOne;
            freqOne = freqTwo;
            freqTwo = temp;
        }

        System.out.println("Freq1: " + freqOne + " Freq2: " + freqTwo);


        RationalApproximation.setAlphNum(freqOne * 100 * 100); // alpha = alph_num / denum,,, = 200 * 100
        RationalApproximation.setDenum(freqTwo * 100 * 100); // = 100 * 401
        RationalApproximation.setDNum(freqOne * 2 * 100); // d = d_num / denum,,, = 401
        int[] fraction = RationalApproximation.mainThing();
        numerator = fraction[0];
        denominator = fraction[1];

        System.out.println("Numerator: " + numerator + " Denominator: " + denominator);

        return numerator + denominator;
    }

    // Credit for Rational Approximation code: https://link.springer.com/chapter/10.1007/978-3-540-72914-3_15
    private static class RationalApproximation {
        private static int alpha_num, d_num, denum;

        private static void setAlphNum(int a) { alpha_num = a; }

        private static void setDenum(int d) { denum = d; }

        private static void setDNum(int d) { d_num = d; }

        private static boolean less(int a,int b,int c,int d) { return (a*d < b*c); }
        private static boolean less_or_equal(int a, int b, int c, int d) { return (a*d <= b*c); }

        // check whether a/b is a valid approximation 
        private static int matches(int a, int b) {
            if (less_or_equal(a,b,alpha_num-d_num,denum)) return 0; 
            if (less(a,b,alpha_num+d_num,denum)) return 1;
            return 0;
        }

        // set initial bounds for the search: 
        static int p_a = 0 ; 
        static int q_a = 1 ; 
        static int p_b = 1 ; 
        static int q_b = 1;

        private static int[] find_exact_solution_left(int p_a, int q_a, int p_b, int q_b) { 
            int k_num = denum * p_b - (alpha_num + d_num) * q_b;
            int k_denum = (alpha_num + d_num) * q_a - denum * p_a;
            int k = (k_num / k_denum) + 1;
            int[] result = {(p_b + k*p_a), (q_b + k*q_a)}; 
            return result;
        }
        private static int[] find_exact_solution_right(int p_a,int q_a,int p_b,int q_b) { 
            int k_num = - denum * p_b + (alpha_num - d_num) * q_b;
            int k_denum = - (alpha_num - d_num) * q_a + denum * p_a;
            int k = (k_num / k_denum) + 1;
            int[] result = {(p_b + k*p_a), (q_b + k*q_a)}; 
            return result;
        }

        private static int[] mainThing() {
            while (true) {
                // Compute the number of steps to the left
                int x_num = denum * p_b - alpha_num * q_b;
                int x_denum = -denum * p_a + alpha_num * q_a;
                int x = (x_num + x_denum - 1) / x_denum; // Equivalent to ceil(x_num / x_denum)

                // Check whether we have a valid approximation
                boolean aa = matches(p_b + x * p_a, q_b + x * q_a) == 1;
                boolean bb = matches(p_b + (x - 1) * p_a, q_b + (x - 1) * q_a) == 1;
                if (aa || bb) {
                    return find_exact_solution_left(p_a, q_a, p_b, q_b);
                }

                // Update the interval
                int new_p_a = p_b + (x - 1) * p_a;
                int new_q_a = q_b + (x - 1) * q_a;
                int new_p_b = p_b + x * p_a;
                int new_q_b = q_b + x * q_a;
                p_a = new_p_a;
                q_a = new_q_a;
                p_b = new_p_b;
                q_b = new_q_b;

                // Compute the number of steps to the right
                x_num = alpha_num * q_b - denum * p_b;
                x_denum = -alpha_num * q_a + denum * p_a;
                x = (x_num + x_denum - 1) / x_denum; // Equivalent to ceil(x_num / x_denum)

                // Check whether we have a valid approximation
                aa = matches(p_b + x * p_a, q_b + x * q_a) == 1;
                bb = matches(p_b + (x - 1) * p_a, q_b + (x - 1) * q_a) == 1;
                if (aa || bb) {
                    return find_exact_solution_right(p_a, q_a, p_b, q_b);
                }

                // Update the interval
                new_p_a = p_b + (x - 1) * p_a;
                new_q_a = q_b + (x - 1) * q_a;
                new_p_b = p_b + x * p_a;
                new_q_b = q_b + x * q_a;
                p_a = new_p_a;
                q_a = new_q_a;
                p_b = new_p_b;
                q_b = new_q_b;
                
                
            }
        }

    }

    /**
     * Creates a file name unique to this pair of player types
     * For example, two random players would result in  "Pair: Random Random.txt"
    **/
    private static String getFileName(Player p1, Player p2) throws Exception {
        String fileName = "Pair_ ";
        fileName+=getPlayerBasedFileID(p1);
        fileName+=getPlayerBasedFileID(p2);
        return fileName + ".txt";
    }

    private static String getPlayerBasedFileID(Player p) throws Exception {
        if (p instanceof RandomPlayer) {
            return "Random";
        }
        else if (p instanceof SimpleReinforcementPlayer) {
            return "Simple";
        }
        else if (p instanceof ChordPlayer) {
            return "Chord";
        }
        else if (p instanceof ScalePlayer) {
            return "Scale";
        }
        else if (p instanceof PredictiveHarmonyPlayer) {
            return "Predictive";
        }
        else if (p instanceof StepwisePlayer) {
            return "Stepwise";
        }
        else if (p instanceof ChordFollowingReinforcementLearning) {
            return "ChordReinforce";
        }
        else if (p instanceof ChordSpecificReinforcementPlayer) {
            return "Measure";
        }
        else if (p instanceof DoubleReinforcementPlayer) {
            return "Double";
        }
        throw new Exception("There's a player type without a file ID");
    }
}
