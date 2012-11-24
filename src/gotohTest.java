import gotoh.GotohAnswer;
import gotoh.GotohProfile;
import gotoh.LocalAligner;

import java.io.IOException;

import loader.Converter;
import sscc.SsccFile;
import sscc.SsccReader;

public class gotohTest {
	public static void main(String[] args) throws IOException {
		GotohProfile prof = new GotohProfile();
		prof.setPrintali(true);
		Converter c = new Converter();

		String[] sequ1 = new String[2];

		sequ1 = "1jcdB00:SSNAKADQASSDAQTANAKADQASNDANAARSDAQAAKDDAARANQRADNA"
				.split(":");

		int[] seq1 = c.convertSeq(sequ1[1]);
		SsccReader reader = new SsccReader("1jy3Q00");
		SsccFile sscc = reader.findAndRead();
		LocalAligner al = new LocalAligner(prof, seq1, sscc, sequ1[0]);
		GotohAnswer ga = new GotohAnswer();
		ga = al.alignPair();
		ga.printAlignment();
	}
}
