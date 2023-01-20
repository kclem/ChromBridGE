ChromBridGE is a utility for detecting and characterizing translocations that arise from genome editing, particularly from CRISPR/Cas cleavage.

Usage:

```
python ChromBridGE.py
options:
  -h, --help            show this help message and exit
  -f FASTQ, --fastq FASTQ
                        Input fastq file
  -a SEQUENCE_A, --sequence_a SEQUENCE_A
                        Input sequence a
  -b SEQUENCE_B, --sequence_b SEQUENCE_B
                        Input sequence b
  --seqA_cut_pos SEQA_CUT_POS
                        Index in sequence a of predicted cut site
  --seqB_cut_pos SEQB_CUT_POS
                        Index in sequence b of predicted cut site
  --match_score MATCH_SCORE
                        Match score for alignment
  --mismatch_score MISMATCH_SCORE
                        Mismatch score for alignment
  --gap_score GAP_SCORE
                        Gap score for alignment
  --jump_score JUMP_SCORE
                        Jump score for alignment
  --cut_pos_incentive_score CUT_POS_INCENTIVE_SCORE
                        Incentive for jumping at a predicted cut site
  -o OUTPUT_FILE, --output_file OUTPUT_FILE
                        Output file to write results
```
