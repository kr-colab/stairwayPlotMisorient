import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.StringTokenizer;

public class Stairway_plot_output_summary_commandline
{
  public static void main(String[] args)
    throws Exception
  {
    if (args.length != 4)
    {
      System.out.println("Usage: java Stairway_plot_output_summary_commandline folder_to_allTheta_files mutation_rate year_per_generation output_file_name ");
      System.out.println("    folder_to_allTheta_files: the folder containing all (and only) the .addTheta files");
      System.out.println("    mutation_rate: assumed mutation rate per site per generation, e.g. 1.2e-8.");
      System.out.println("    year_per_generation: assumed generation time (in years), e.g. 24.");
      System.out.println("    output_file_name: name for the output summary file");
      System.exit(1);
    }
    String dir = correctDirectory(args[0]);
    File folder = new File(dir);
    File[] listOfFiles = folder.listFiles();
    String[] infiles = folder.list();
    Arrays.sort(infiles);
    double mu = Double.parseDouble(args[1]);
    double year_per_generation = Double.parseDouble(args[2]);
    String outfile = args[3];
    PrintWriter out = new PrintWriter(new FileWriter(outfile), true);
    
    BufferedReader in = new BufferedReader(new FileReader(dir + infiles[0]));
    String line = in.readLine();
    out.println(line);
    
    StringTokenizer t = new StringTokenizer(line, "\t");
    int nt = t.countTokens();
    
    String[] mspara = new String[nt];
    for (int i = 0; i < nt; i++) {
      mspara[i] = t.nextToken();
    }
    int nseq = Integer.parseInt(mspara[1]);
    int nrep = infiles.length;
    long L = Long.parseLong(mspara[2]);
    
    System.out.println("nrep=" + nrep);
    out.println("nrep=" + nrep);
    double[][] allest = new double[nseq - 1][nrep];
    out.println("final fitness:");
    for (int ii = 0; ii < nrep - 1; ii++) {
      out.print(infiles[ii] + "\t");
    }
    out.println(infiles[(nrep - 1)]);
    for (int ii = 0; ii < nrep; ii++)
    {
      in = new BufferedReader(new FileReader(dir + infiles[ii]));
      in.readLine();
      in.readLine();
      System.out.println(ii + 1 + "\t" + infiles[ii]);
      line = in.readLine();
      double fitness = 0.0D;
      while (line.indexOf("final model:") == -1) {
        line = in.readLine();
      }
      t = new StringTokenizer(line, "\t");
      t.nextToken();
      fitness = Double.parseDouble(t.nextToken());
      if (ii == nrep - 1) {
        out.print(fitness);
      } else {
        out.print(fitness + "\t");
      }
      line = in.readLine();
      t = new StringTokenizer(line);
      int dim = t.countTokens();
      int[][] group = new int[dim][];
      for (int i = 0; i < dim; i++)
      {
        String tmp = t.nextToken();
        StringTokenizer t2 = new StringTokenizer(tmp, ",");
        int nxi = t2.countTokens();
        group[i] = new int[nxi];
        for (int j = 0; j < nxi; j++) {
          group[i][j] = Integer.parseInt(t2.nextToken());
        }
      }
      line = in.readLine();
      t = new StringTokenizer(line);
      for (int i = 0; i < dim; i++)
      {
        double theta = Double.parseDouble(t.nextToken());
        for (int j = 0; j < group[i].length; j++) {
          allest[(group[i][j] - 2)][ii] = theta;
        }
      }
      in.close();
    }
    out.println();
    int dim = nseq - 1;
    double[] time = new double[dim];
    for (int i = 0; i < dim; i++)
    {
      Arrays.sort(allest[i]);
      time[i] = (allest[i][(nrep / 2)] / L / (i + 2) / (i + 1));
    }
    for (int i = dim - 2; i >= 0; i--) {
      time[i] += time[(i + 1)];
    }
    out.println("mutation_per_site\ttheta\ttheta_per_site_median\ttheta_per_site_2.5%\ttheta_per_site_97.5%\tyear\tNe_median\tNe_2.5%\tNe_97.5%");
    for (int i = dim - 1; i >= 0; i--) {
      if (i == dim - 1)
      {
        out.print("0\t" + (i + 2) + "\t" + allest[i][(nrep / 2)] / L + "\t" + allest[i][((int)(nrep * 0.025D))] / L + "\t" + allest[i][((int)(nrep * 0.975D))] / L + "\t" + 0.0D / mu * year_per_generation + "\t" + allest[i][(nrep / 2)] / L / mu / 4.0D + "\t" + allest[i][((int)(nrep * 0.025D))] / L / mu / 4.0D + "\t" + allest[i][((int)(nrep * 0.975D))] / L / mu / 4.0D);
        out.println();
        out.print(time[i] + "\t" + (i + 2) + "\t" + allest[i][(nrep / 2)] / L + "\t" + allest[i][((int)(nrep * 0.025D))] / L + "\t" + allest[i][((int)(nrep * 0.975D))] / L + "\t" + time[i] / mu * year_per_generation + "\t" + allest[i][(nrep / 2)] / L / mu / 4.0D + "\t" + allest[i][((int)(nrep * 0.025D))] / L / mu / 4.0D + "\t" + allest[i][((int)(nrep * 0.975D))] / L / mu / 4.0D);
        out.println();
      }
      else
      {
        out.print(time[(i + 1)] + "\t" + (i + 2) + "\t" + allest[i][(nrep / 2)] / L + "\t" + allest[i][((int)(nrep * 0.025D))] / L + "\t" + allest[i][((int)(nrep * 0.975D))] / L + "\t" + time[(i + 1)] / mu * year_per_generation + "\t" + allest[i][(nrep / 2)] / L / mu / 4.0D + "\t" + allest[i][((int)(nrep * 0.025D))] / L / mu / 4.0D + "\t" + allest[i][((int)(nrep * 0.975D))] / L / mu / 4.0D);
        out.println();
        out.print(time[i] + "\t" + (i + 2) + "\t" + allest[i][(nrep / 2)] / L + "\t" + allest[i][((int)(nrep * 0.025D))] / L + "\t" + allest[i][((int)(nrep * 0.975D))] / L + "\t" + time[i] / mu * year_per_generation + "\t" + allest[i][(nrep / 2)] / L / mu / 4.0D + "\t" + allest[i][((int)(nrep * 0.025D))] / L / mu / 4.0D + "\t" + allest[i][((int)(nrep * 0.975D))] / L / mu / 4.0D);
        out.println();
      }
    }
    out.close();
  }
  
  private static String correctDirectory(String dir)
  {
    StringTokenizer t = new StringTokenizer(dir, "/\\\"'");
    int nt = t.countTokens();
    String sep = File.separator;
    String dir2 = "";
    for (int i = 0; i < nt; i++) {
      dir2 = dir2 + t.nextToken() + sep;
    }
    if (((dir.charAt(0) == '"') || (dir.charAt(0) == '\'')) && ((dir.charAt(1) == '/') || (dir.charAt(1) == '\\'))) {
      dir2 = sep + dir2;
    } else if ((dir.charAt(0) == '/') || (dir.charAt(0) == '\\')) {
      dir2 = sep + dir2;
    }
    return dir2;
  }
}
