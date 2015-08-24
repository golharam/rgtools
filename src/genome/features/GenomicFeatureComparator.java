package genome.features;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

/* 
 * Every genomic feature can be compared using chromosome and start location.
 */
public class GenomicFeatureComparator implements Comparator<GenomicFeature> {
    /**
     * Create a new list that is ordered by the starting index of the features' locations. 
     *
     * @return An ordered list of feature based on chromosome and start positions based on the positive strand.
     */
    public static List<GenomicFeature> sortByGenomicStart(List<GenomicFeature> features) {
    	GenomicFeature array[] = features.toArray(new GenomicFeature[1]);

        Arrays.sort(array, new GenomicFeatureComparator());

        return new ArrayList<GenomicFeature>(Arrays.asList(array));
    }

    public int compare(GenomicFeature a, GenomicFeature b) {
		char a_chr = a.getChromosome().charAt(3);
		char b_chr = b.getChromosome().charAt(3);
		if (a_chr < b_chr) {
			return -1;
		}
		if (a_chr > b_chr) {
			return 1;
		}
		if (a_chr == b_chr) {
			return a.getStart() - b.getStart();
		} else {
			throw new IndexOutOfBoundsException(
					"Something went wrong comparing two features:\n"
							+ a.toString() + "\n" + b.toString());
		}
	}

    
//    public int compare(GFFFeature a, GFFFeature b) {
 //   	return compare((GenomicFeatureI)a, (GenomicFeatureI)b);
  //  }
}

