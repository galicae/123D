package loader;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

/**
 * this class handles all operations that have to do with reading files; it
 * should take care of reading the pairFile and retrieve the aforementioned
 * pairs from a sequence library file.
 * 
 * @author nikos papadopoulos
 * 
 */
public class JLoader {
	private String pairFile;
	private String seqLibFile;
	public String[][] pairs;
	public String[] sequences;
	private int pairLength;
	private boolean ordered = true;

	public JLoader(String pairFile, String seqLibFile) {
		this.pairFile = pairFile;
		this.seqLibFile = seqLibFile;
	}

	/**
	 * this function parses a .pairs file and stores the first two columns of
	 * each row, the IDs of the sequences to be aligned
	 * 
	 * @param pairFile
	 *            a file containing at least two columns per row, containing the
	 *            IDs of sequences that are to be aligned with one another
	 * @return a String[][] with each array being a pair of IDs that need to be
	 *         aligned
	 */
	public String[][] loadPairFile() {
		try {
			FileInputStream fstream = new FileInputStream(pairFile);
			DataInputStream in = new DataInputStream(fstream);
			BufferedReader br = new BufferedReader(new InputStreamReader(in));
			String strLine;
			pairLength = countLines(pairFile);
			pairs = new String[pairLength][2];
			int i = 0;
			while ((strLine = br.readLine()) != null) {
				pairs[i][0] = strLine.split(" ")[0];
				pairs[i][1] = strLine.split(" ")[1];
				i++;
			}
			in.close();
			return pairs;
		} catch (Exception e) {// Catch exception if any
			System.err.println("Error: The pairfile has this to say: " + e.getMessage());
		}
		return null;
	}

	/**
	 * *really* fast search for \n characters - maybe expand for whole
	 * expressions?
	 * 
	 * @param filename
	 *            the name of the file to count lines into
	 * @return the number of lines of the text file
	 * @throws IOException
	 */
	public int countLines(String filename) throws IOException {
		InputStream is = new BufferedInputStream(new FileInputStream(filename));
		try {
			byte[] c = new byte[1024];
			int count = 0;
			int readChars = 0;
			boolean empty = true;
			while ((readChars = is.read(c)) != -1) {
				empty = false;
				for (int i = 0; i < readChars; ++i) {
					if (c[i] == '\n')
						++count;
				}
			}
			return (count == 0 && !empty) ? 1 : count;
		} finally {
			is.close();
		}
	}

	/**
	 * this function searches for a given query string in a list of strings. The
	 * assumption is that the query string is the ID of a sequence in the
	 * seqLibFile, which is why the search doesn't match the whole string but
	 * the first 7 characters.
	 * 
	 * @param query
	 *            the ID of the sequence we are looking for
	 * @param seqList
	 *            an array where the sequences from the pair file are listed in
	 *            alphabetical order
	 * @return the sequence with {@code query} as its ID, or null if the query
	 *         cannot be satisfied
	 */
	public String retrieveSeqBin(String query, String[] seqList) {
		/*
		 * binary search. Seemed to be an average of 0.02 secs better than
		 * manual search for queries near the middle of the test file. Since
		 * this procedure will be used twice every alignment for a large total
		 * of sequences , I'm guessing the sacrifice in memory is worth the time
		 * I'll be winning.
		 * 
		 * it remains a very interesting question what would happen with a
		 * really big seqLib file; will a computer have enough memory for all
		 * things to be stored in RAM? Remember I am also saving the pair list
		 * on memory.
		 * 
		 * I think the computer will easily handle this assignment, but I should
		 * probably look for a larger seqlib and pairfile.
		 */
		int min = 0;
		int max = seqList.length;
		while (min < max) {
			int mid = (min + max) / 2;
			assert (mid < max);
			// ah, the compare method for strings... so nice to have libraries!
			int comp = query.compareToIgnoreCase(seqList[mid]);
			if (comp < 0)
				max = mid;
			else if (comp == 0)
				return seqList[mid];
			else
				min = mid + 1;
		}
		if ((max == min) && (seqList[min].startsWith(query))) {
//			System.out.println("Query found at position " + min + ":"
//					+ seqList[min]);
			return seqList[min];
		} else
			System.out.println("Query " + query + " not found.");
		return null;
	}

	/**
	 * this function searches for a given query string in a list of strings. The
	 * assumption is that the query string is the ID of a sequence in the
	 * seqLibFile, which is why the search doesn't match the whole string but
	 * the first 7 characters.
	 * 
	 * @param query
	 *            the ID of the sequence we are looking for
	 * @param seqList
	 *            a list with all sequences contained in the seqLibFile. Needs
	 *            not be in alphabetical order
	 * @return the sequence with {@code query} as its ID, or null if the query
	 *         cannot be satisfied
	 */
	public String retrieveSeqMan(String query, String[] seqList) {
		int i;
		for (i = 0; i < seqList.length; i++) {
			if (!seqList[i].startsWith(query))
				continue;
			else
				break;
		}
		if (i < seqList.length)
			return seqList[i];
		else
			System.out.println("Query " + query + " not found.");
		return null;
	}

	/**
	 * this function reads a seqlib-File and stores it as an array
	 * 
	 * @param file
	 *            the seqlib-File to load
	 * @return an array with sequences at the cells
	 * @throws IOException
	 */
	public String[] loadSeqLibFile() throws IOException {
		int lines = this.countLines(seqLibFile);
		sequences = new String[lines];
		FileInputStream fstream = new FileInputStream(seqLibFile);
		DataInputStream in = new DataInputStream(fstream);
		BufferedReader br = new BufferedReader(new InputStreamReader(in));
		String strLine = "";
		String prev = "";
		int i = 0;
		// save each sequence at a large list
		while ((strLine = br.readLine()) != null) {
//			if (ordered)
//				if (strLine.subSequence(0, 6).toString().compareTo(prev) < 0)
//					ordered = false;
			sequences[i] = strLine;
			prev = strLine;
			i++;
		}
		// remember to close the stream!
		in.close();
		return sequences;
	}

	/**
	 * this function retrieves a pair of sequences from the sequence array and
	 * returns them for alignment
	 * 
	 * @param seq1
	 *            the first sequence
	 * @param seq2
	 *            the second sequence
	 * @param seqList
	 *            the array of sequences
	 * @return a String[2] array with sequences 1 and 2
	 */
	public String[] retrieveSeqPair(String seq1, String seq2, String[] seqList) {
		// now is the time to see if the list was ordered or not!
		seqList = new String[2];
		if (ordered) {
			seqList[0] = retrieveSeqBin(seq1, seqList);
			seqList[1] = retrieveSeqBin(seq2, seqList);
		} else {
			seqList[0] = retrieveSeqMan(seq1, seqList);
			seqList[1] = retrieveSeqMan(seq2, seqList);
		}
		return seqList;
	}
	
	
	public int getPairLength() {
		return pairLength;
	}

}
