import java.io.IOException;

import PSS.PssReader;

public class pssTest {
	public static void main(String[] args) throws IOException {
		String dir = "sscc/";
		PssReader reader = new PssReader(dir);
		reader.calcSecStructPref();
		System.out.println();
	}
}
