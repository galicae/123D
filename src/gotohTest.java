import gotoh.FreeshiftAligner;
import gotoh.GotohAnswer;
import gotoh.GotohProfile;

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

		sequ1 = "1wq2B00:MEEAKQKVVDFLNSKSKSKFYFNDFTDLFPDMKQREVKKILTALVNDEVLEYWSSGSTTMYGLKG"
				.split(":");

		int[] seq1 = c.convertSeq(sequ1[1]);
		SsccReader reader = new SsccReader("1j2xA00");
		SsccFile sscc = reader.findAndRead();
		FreeshiftAligner al = new FreeshiftAligner(prof, seq1, sscc, "1wq2B00");
		GotohAnswer ga = new GotohAnswer();
		ga = al.alignPair();
		ga.printAlignment();
	}

}
