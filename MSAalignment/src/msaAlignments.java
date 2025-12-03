import java.util.*;
import java.io.*;
public class msaAlignments {

    //this contains two different ways to align sequences, the normal way we learned in class
    //and my understanding/take of the A* algorithm (from a paper here: https://cdn.aaai.org/AAAI/2002/AAAI02-111.pdf)
    //i then have ways to check the differences in memory and time usage, along with if it has the optimal alignment so all this can be used in analysis

    private static final int MATCH = 2; //matches are +2
    private static final int MISMATCH = -1; //mismatches are -1
    private static final int GAP = -2; //gaps are -2

    // classes to hold data, this one is to hold the sequence data from the files
    static class Sequence {
        String id; //id of the sequence, uniprot id!!
        String seq; //dna seq

        // actual constructor for the sequence
        Sequence(String id, String seq) {
            this.id = id;
            this.seq = seq.toUpperCase();//uppercase so everything is the same no matter what
        }

        //lengh method yes
        int length() {
            return seq.length();
        }
    }

    // will hold the results of the alignments
    static class AlignmentResult {
        String aseq1;//first sequence with gaps
        String aseq2;//second sequence with gaps
        int score;//total alignment score
        long time;//time in nanoseconds
        long memory; //memory in bytes so it can be compared
        int cellsComputed; //how many dp cells calculated
        int nodesExpanded; //how many nodes look at (only used for my take on A*)

        // another constructor
        AlignmentResult(String a1, String a2, int score) {
            this.aseq1 = a1;
            this.aseq2 = a2;
            this.score = score;
        }
    }

    // has comparison of the results between the a=two algorithms
    static class ComparisonResult {
        String seq1Id, seq2Id; //sequence ids
        int len1, len2;//sequence lengths

        //the standard algo DP results
        int dpScore;
        long dpTime;
        long dpMemory;
        int dpCells;

        //A*-pruned results
        int astarScore;
        long astarTime;
        long astarMemory;
        int astarCells;
        int astarNodes;

        boolean scoresMatch; //boolean asking if scores match (because i'm also testing if they have optimal alignments)
        double speedRatio;//a*Time / dpTime
        double cellsRatio; // a*Cells / dpCells

        // another constructor
        ComparisonResult(String id1, String id2, int l1, int l2) {
            this.seq1Id = id1;
            this.seq2Id = id2;
            this.len1 = l1;
            this.len2 = l2;
        }

        //gets the ratios between the two after alignments are done
        void calculateRatios() {
            scoresMatch = (dpScore == astarScore); //gets the boolean val
            speedRatio = (double) astarTime / dpTime;
            cellsRatio = (double) astarCells / dpCells;
        }
    }

    public static class StandardDP { //this is the standard dp tabke program that we learned

        public static AlignmentResult align(String s1, String s2) { //
            Runtime runtime = Runtime.getRuntime(); //start runtime count
            System.gc();//try to clean system memory so i'm only calculating this algo's memory
            long startMemory = runtime.totalMemory() - runtime.freeMemory(); //try to get memory
            long startTime = System.nanoTime();//current time (ns)
            int m = s1.length(); //length of first seq
            int n = s2.length(); //length of second seq

            //dp table using the normal bellman equation for an alignment
            int[][] dp = new int[m + 1][n + 1];
            int cellsComputed = 0;//keep counting how many cells that are getting looked at

            //make the first row and column (always all gaps)
            for (int i = 0; i <= m; i++) {
                dp[i][0] = i * GAP;
                cellsComputed++; //count the cells (even if they're base cases i'm still counting them)
            }
            for (int j = 0; j <= n; j++) {
                dp[0][j] = j * GAP;
                cellsComputed++; //count the cells (even if they're base cases i'm still counting them)
            }

            //fill the whole dp table
            for (int i = 1; i <= m; i++) {
                for (int j = 1; j <= n; j++) {
                    // match or mismatch
                    int matchScore = dp[i - 1][j - 1] + (s1.charAt(i - 1) == s2.charAt(j - 1) ? MATCH : MISMATCH);
                    // gap in the second string
                    int deleteScore = dp[i - 1][j] + GAP;
                    // gap in the first string
                    int insertScore = dp[i][j - 1] + GAP;
                    dp[i][j] = Math.max(matchScore, Math.max(deleteScore, insertScore)); //dp table always takes best score
                    cellsComputed++; //keep adding on cells computed
                }
            }
            // backtrack through the table (always starting form bottom right)
            StringBuilder aligned1 = new StringBuilder();
            StringBuilder aligned2 = new StringBuilder();
            int i = m, j = n;
            while (i > 0 || j > 0) {
                if (i > 0 && j > 0 && dp[i][j] == dp[i - 1][j - 1] + (s1.charAt(i - 1) == s2.charAt(j - 1) ? MATCH : MISMATCH)) { //match OR a mismatch (diagonal)
                    aligned1.append(s1.charAt(i - 1));
                    aligned2.append(s2.charAt(j - 1));
                    i--;
                    j--;
                } else if (i > 0 && dp[i][j] == dp[i - 1][j] + GAP) { //gap in second string
                    aligned1.append(s1.charAt(i - 1));
                    aligned2.append('-');
                    i--;
                } else { //gap in this first string
                    aligned1.append('-');
                    aligned2.append(s2.charAt(j - 1));
                    j--;
                }
            }

            //both aligned backwards, so they must be reversed
            //we have two because an alignment is sequences on top of one another (which i'll print later)
            aligned1.reverse();
            aligned2.reverse();

            //get final time and memory for the normal dp table
            long endTime = System.nanoTime();
            long endMemory = runtime.totalMemory() - runtime.freeMemory();

            //make the result
            AlignmentResult result = new AlignmentResult(aligned1.toString(), aligned2.toString(), dp[m][n]);  //opt is bottom right as always
            result.time = endTime - startTime;
            result.memory = Math.max(0, endMemory - startMemory);//cant have negatives
            result.cellsComputed = cellsComputed;
            return result;
        }
    }

    public static class AStarPrunedDP {
        //a* heuristic pruned way of searching, it only computed the dp table cells it NEEDS to
        //the huge thing about this one is that because its a greedy algorithm and might not always get the optimal alignment

        static class Node implements Comparable<Node> {
            int i, j;//position in DP table (row, column)
            int g;// actual score so far (from DP table)
            int h; //heurtsic score estimate (to see what we need to calc)
            Node parent; // for backtracking

            //constructor
            Node(int i, int j, int g, int h, Node parent) { //just with the constructor
                this.i = i;
                this.j = j;
                this.g = g;
                this.h = h;
                this.parent = parent;
            }

            //total score f = g + h (A* formula)
            int f() {
                return g + h;
            }

            // compare nodes to see priority, i want nodes with HIGHER f-score
            @Override
            public int compareTo(Node other) {
                return Integer.compare(other.f(), this.f());
            }
        }

        //heuristic function where it estimates best possible score
        private static int heuristic(String s1, String s2, int i, int j) {
            int remaining1 = s1.length() - i;//chars left in string 1
            int remaining2 = s2.length() - j; //chars left in string 2

            //best case wjere all remaining characters match
            int maxMatches = Math.min(remaining1, remaining2);
            int unavoidableGaps = Math.abs(remaining1 - remaining2);
            return maxMatches * MATCH + unavoidableGaps * GAP;
        }

        public static AlignmentResult align(String s1, String s2) {
            Runtime runtime = Runtime.getRuntime(); //start runtime
            System.gc();
            long startMemory = runtime.totalMemory() - runtime.freeMemory(); //start memory
            long startTime = System.nanoTime();
            int m = s1.length(); //get string length
            int n = s2.length();//get string length

            // DP table but we'll compute cells only when needed
            int[][] dp = new int[m + 1][n + 1];
            boolean[][] computed = new boolean[m + 1][n + 1];  //track computed cells
            int cellsComputed = 0;//count cells that we've been through

            //initialize base cases
            for (int i = 0; i <= m; i++) {
                dp[i][0] = i * GAP;
                computed[i][0] = true;
                cellsComputed++;
            }
            for (int j = 0; j <= n; j++) {
                dp[0][j] = j * GAP;
                computed[0][j] = true;
                cellsComputed++;
            }

            //a* search setup
            PriorityQueue<Node> openSet = new PriorityQueue<>();//nodes to look through
            Set<String> closedSet = new HashSet<>();//nodes we've been to
            int nodesExpanded = 0;//how many nodes we've seen

            //start at top left node
            openSet.add(new Node(0, 0, dp[0][0], heuristic(s1, s2, 0, 0), null));
            Node goalNode = null;//set this when the end is reached

            //the main a* loop
            while (!openSet.isEmpty()) { //get the node with the best score
                Node current = openSet.poll();
                nodesExpanded++;

                if (current.i == m && current.j == n) { //check if the goal has been reached (bottom right)
                    goalNode = current;
                    break;
                }

                String key = current.i + "," + current.j;
                if (closedSet.contains(key)) {
                    continue;
                }
                closedSet.add(key); //mark as visited

                //make sure current cell has been computed
                if (!computed[current.i][current.j]) {
                    if (current.i > 0 && current.j > 0) {//
                        int match = dp[current.i - 1][current.j - 1] +
                                (s1.charAt(current.i - 1) == s2.charAt(current.j - 1) ? MATCH : MISMATCH);
                        int delete = dp[current.i - 1][current.j] + GAP;
                        int insert = dp[current.i][current.j - 1] + GAP;
                        dp[current.i][current.j] = Math.max(match, Math.max(delete, insert));
                        computed[current.i][current.j] = true;
                        cellsComputed++; //add to the num of computed cells
                    }
                }

                if (current.i < m) { //move down (gap in seq 2)
                    int nextI = current.i + 1;
                    if (!computed[nextI][current.j]) { //compute cell if we haven;t alreasdy
                        dp[nextI][current.j] = dp[nextI - 1][current.j] + GAP;
                        computed[nextI][current.j] = true;
                        cellsComputed++;
                    }

                    //add to the open set
                    openSet.add(new Node(nextI, current.j, dp[nextI][current.j], heuristic(s1, s2, nextI, current.j), current));
                }

                if (current.j < n) { //move ot the right (gap in seq 1)
                    int nextJ = current.j + 1;
                    if (!computed[current.i][nextJ]) {
                        dp[current.i][nextJ] = dp[current.i][nextJ - 1] + GAP;
                        computed[current.i][nextJ] = true;
                        cellsComputed++;
                    }

                    //add to open set
                    openSet.add(new Node(current.i, nextJ, dp[current.i][nextJ], heuristic(s1, s2, current.i, nextJ), current));
                }

                if (current.i < m && current.j < n) { //move diagonally (so match or mismatch)
                    int nextI = current.i + 1;
                    int nextJ = current.j + 1;

                    if (!computed[nextI][nextJ]) { //need to compute the diagonal cell
                        int match = dp[current.i][current.j] +
                                (s1.charAt(current.i) == s2.charAt(current.j) ? MATCH : MISMATCH);
                        int delete = dp[nextI][current.j] + GAP;
                        int insert = dp[current.i][nextJ] + GAP;
                        dp[nextI][nextJ] = Math.max(match, Math.max(delete, insert));
                        computed[nextI][nextJ] = true;
                        cellsComputed++; //add to num of cells computed
                    }

                    //once again, add to the open set
                    openSet.add(new Node(nextI, nextJ, dp[nextI][nextJ], heuristic(s1, s2, nextI, nextJ), current));
                }
            }

            //make the alignments with the a* algo
            StringBuilder aligned1 = new StringBuilder();
            StringBuilder aligned2 = new StringBuilder();

            if (goalNode != null) { //traceback through the pointers
                Node current = goalNode;
                while (current.parent != null) { //figure out what move was made
                    int di = current.i - current.parent.i;//chage in te row
                    int dj = current.j - current.parent.j; //change in the column

                    if (di == 1 && dj == 1) { //diagonal move
                        aligned1.append(s1.charAt(current.i - 1));
                        aligned2.append(s2.charAt(current.j - 1));
                    } else if (di == 1 && dj == 0) { //down move
                        aligned1.append(s1.charAt(current.i - 1));
                        aligned2.append('-');
                    } else { //right move
                        aligned1.append('-');
                        aligned2.append(s2.charAt(current.j - 1));
                    }
                    current = current.parent;
                }
                //once again it was built backwards so it has to be reversed
                aligned1.reverse();
                aligned2.reverse();
            }

            long endTime = System.nanoTime(); //time
            long endMemory = runtime.totalMemory() - runtime.freeMemory(); //runtime
            AlignmentResult result = new AlignmentResult( //make new alignment
                    aligned1.toString(),
                    aligned2.toString(),
                    goalNode != null ? goalNode.g : 0);
            result.time = endTime - startTime;
            result.memory = Math.max(0, endMemory - startMemory);
            result.cellsComputed = cellsComputed;
            result.nodesExpanded = nodesExpanded;
            return result;
        }
    }

    public static List<Sequence> parseFasta(String filename) throws IOException { //reads sequences from a fasta file
        List<Sequence> sequences = new ArrayList<>();
        BufferedReader reader = new BufferedReader(new FileReader(filename));
        String line;
        String currentId = null;
        StringBuilder currentSequence = new StringBuilder();
        while ((line = reader.readLine()) != null) {
            line = line.trim();//remove the end spaces
            if (line.isEmpty()) continue; //skip the empty lines (theres sometimes breaks between seqs)
            if (line.startsWith(">")) { //header always starts with > so look for this
                if (currentId != null) {
                    sequences.add(new Sequence(currentId, currentSequence.toString()));
                }
                String[] parts = line.substring(1).split("\\s+", 2);//split on white space
                currentId = parts[0];//first part is the id
                currentSequence = new StringBuilder();//start the new sewquence
            } else { //sequence line
                currentSequence.append(line.toUpperCase().replaceAll("\\s+", ""));
            }
        }
        if (currentId != null) { //;ast sequence
            sequences.add(new Sequence(currentId, currentSequence.toString()));
        }
        reader.close();
        return sequences;
    }

    //here's where the comparisons between the two starts
    public static void compareAlgorithms(List<Sequence> sequences, String outputFile) throws IOException {
        PrintWriter writer = new PrintWriter(new FileWriter(outputFile));
        //write what we're looking at
        writer.printf("%-10s %-10s %6s %6s %10s %10s %10s %10s %10s %10s %10s %8s %8s\n", "Seq1", "Seq2", "Len1", "Len2", "DP Time", "A* Time", "DP Mem", "A* Mem", "DP Cells", "A* Cells", "A* Nodes", "Speed", "Cells");
        List<ComparisonResult> allResults = new ArrayList<>(); //comparison results
        int totalPairs = sequences.size() * (sequences.size() - 1) / 2;
        int currentPair = 0;
        //compare each pair of sequences
        for (int i = 0; i < sequences.size(); i++) {
            for (int j = i + 1; j < sequences.size(); j++) {
                currentPair++;
                Sequence seq1 = sequences.get(i);
                Sequence seq2 = sequences.get(j);
                System.out.printf("[%d/%d] Comparing %s vs %s...\n",
                        currentPair, totalPairs, seq1.id, seq2.id);
                //make the comparison object
                ComparisonResult comp = new ComparisonResult(
                        seq1.id, seq2.id, seq1.length(), seq2.length());
                //run a standard dp table
                System.out.print("standard dp");
                AlignmentResult dpResult = StandardDP.align(seq1.seq, seq2.seq);
                comp.dpScore = dpResult.score;
                comp.dpTime = dpResult.time;
                comp.dpMemory = dpResult.memory;
                comp.dpCells = dpResult.cellsComputed;
                System.out.printf("%,d μs, %,d cells\n", dpResult.time / 1000, dpResult.cellsComputed);
                //run the A* pruned DP table
                System.out.print("pruned table");
                AlignmentResult astarResult = AStarPrunedDP.align(seq1.seq, seq2.seq);
                comp.astarScore = astarResult.score;
                comp.astarTime = astarResult.time;
                comp.astarMemory = astarResult.memory;
                comp.astarCells = astarResult.cellsComputed;
                comp.astarNodes = astarResult.nodesExpanded;
                System.out.printf("%,d μs, %,d cells, %,d nodes\n", astarResult.time / 1000, astarResult.cellsComputed, astarResult.nodesExpanded);
                //calc the ratios
                comp.calculateRatios();
                allResults.add(comp);
                //write it to the files
                writer.printf("%-10s %-10s %6d %6d %,10d %,10d %,10d %,10d %,10d %,10d %,10d %8.2f %8.2f\n", seq1.id, seq2.id, seq1.length(), seq2.length(), dpResult.time / 1000, astarResult.time / 1000, dpResult.memory / 1024, astarResult.memory / 1024, dpResult.cellsComputed, astarResult.cellsComputed,
                        astarResult.nodesExpanded, comp.speedRatio, comp.cellsRatio);
                // chekc if scores match
                if (!comp.scoresMatch) {
                    System.out.println("scores don't match");
                    writer.println("scorea dont match");
                } else {
                    System.out.printf("scores match)\n", dpResult.score);
                }
                System.out.println();
            }
        }
        //calc summary stats
        long totalDpTime = 0, totalAstarTime = 0;
        long totalDpMemory = 0, totalAstarMemory = 0;
        int totalDpCells = 0, totalAstarCells = 0;
        int correctAlignments = 0;
        for (ComparisonResult r : allResults) {
            totalDpTime += r.dpTime;
            totalAstarTime += r.astarTime;
            totalDpMemory += r.dpMemory;
            totalAstarMemory += r.astarMemory;
            totalDpCells += r.dpCells;
            totalAstarCells += r.astarCells;
            if (r.scoresMatch) correctAlignments++;
        }
        double avgDpTime = totalDpTime / (double) allResults.size() / 1000; // Convert to μs
        double avgAstarTime = totalAstarTime / (double) allResults.size() / 1000;
        double avgDpMemory = totalDpMemory / (double) allResults.size() / 1024; // Convert to KB
        double avgAstarMemory = totalAstarMemory / (double) allResults.size() / 1024;
        double avgDpCells = totalDpCells / (double) allResults.size();
        double avgAstarCells = totalAstarCells / (double) allResults.size();

        writer.printf("Total comparisons: %d\n", allResults.size());
        writer.printf("Correct alignments: %d/%d (%.1f%%)\n", correctAlignments, allResults.size(), (100.0 * correctAlignments) / allResults.size());
        writer.println();
        writer.printf("Average Standard DP time: %,.1f μs\n", avgDpTime);
        writer.printf("Average A*-Pruned time:   %,.1f μs\n", avgAstarTime);
        writer.printf("Average speed ratio (A*/DP): %.2fx\n", avgAstarTime / avgDpTime);
        writer.println();
        writer.printf("Average Standard DP memory: %,.1f KB\n", avgDpMemory);
        writer.printf("Average A*-Pruned memory:   %,.1f KB\n", avgAstarMemory);
        writer.printf("Average memory ratio (A*/DP): %.2fx\n", avgAstarMemory / avgDpMemory);
        writer.println();
        writer.printf("Average Standard DP cells computed: %,.1f\n", avgDpCells);
        writer.printf("Average A*-Pruned cells computed:   %,.1f\n", avgAstarCells);
        writer.printf("Average cells ratio (A*/DP): %.2fx\n", avgAstarCells / avgDpCells);
        writer.printf("Cell reduction: %.1f%%\n", 100 * (1 - avgAstarCells / avgDpCells));
        //best and worst cases
        ComparisonResult fastestAstar = allResults.get(0);
        ComparisonResult slowestAstar = allResults.get(0);
        ComparisonResult bestPruning = allResults.get(0);//lowest cell ratio
        for (ComparisonResult r : allResults) {
            if (r.speedRatio < fastestAstar.speedRatio) fastestAstar = r;
            if (r.speedRatio > slowestAstar.speedRatio) slowestAstar = r;
            if (r.cellsRatio < bestPruning.cellsRatio) bestPruning = r;
        }
        writer.printf("Fastest A*: %s vs %s (%.2fx faster than DP)\n",
                fastestAstar.seq1Id, fastestAstar.seq2Id, 1 / fastestAstar.speedRatio);
        writer.printf("Slowest A*: %s vs %s (%.2fx slower than DP)\n",
                slowestAstar.seq1Id, slowestAstar.seq2Id, slowestAstar.speedRatio);
        writer.printf("Best pruning: %s vs %s (only %.1f%% of DP cells computed)\n",
                bestPruning.seq1Id, bestPruning.seq2Id, 100 * bestPruning.cellsRatio);
        writer.close();
        System.out.println("Results saved to: " + outputFile);
    }

    //where program for comparisons starts
    public static void main(String[] args) {
        System.out.println("Standard DP vs A*-Pruned DP");
        //command line args
        if (args.length == 0) {
            return; }
        String fastaFile = args[0];
        String outputFile = args.length > 1 ? args[1] : "results.txt";
        try {
                //load sequences
                System.out.println("loading from: " + fastaFile);
                List<Sequence> sequences = parseFasta(fastaFile);
                if (sequences.size() < 2) {
                    System.out.println("need at least 2 seqs");
                    return;
                }
                System.out.println("loaded " + sequences.size() + " sequences:"); //to check if you're getting all of themm
                //goes through each seq in list
                for (Sequence s : sequences) {
                    System.out.printf("  %s (length: %d)\n", s.id, s.length());
                }
                //run comparison
                compareAlgorithms(sequences, outputFile);

            } catch (IOException e) {
                System.out.println(e.getMessage());
            }
        }
    }