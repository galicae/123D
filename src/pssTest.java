import java.io.IOException;

import resources.AminoAcids;

import PSS.PssReader;

public class pssTest {
	public static void main(String[] args) throws IOException {
		String dir = "sscc/";
		PssReader reader = new PssReader(dir);
		double[][] result = reader.calcPSSScore();
		for(int i = 0; i < 20; i++) {
			System.out.print(AminoAcids.reverse[i] + "\t");
		}
		System.out.println();
		for (int i = 0; i < 3; i++) {
			for(int j = 0; j < 20; j++) {
				System.out.print((int)(result[j][i]*100) + "\t");
			}
			System.out.println();
		}
//		double logodds = 1.414;
//		System.out.println(Math.log(logodds));
	}
}
