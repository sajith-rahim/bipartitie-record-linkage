*compute distortions*

 * the distortion distribution of x' | y' for field l

for every category value in (M_l), compute it's distortion (string similarity) with every other value in M_l and normalize by sum.

`
{(0, 'SILKE'): defaultdict(<class 'float'>, {'SILKE': 0.04685408299866133, 'PETRA': 0.0, 'RENATE': 0.00780901383311022, 'SUSANNE': 0.01338688085676038,..},
(0, 'PETRA'): defaultdict(<class 'float'>, {'SILKE': 0.0, 'PETRA': 0.04214258240379955,....}
 }
 `
 
 `{
       (FIELD_INDEX, FIELD_VALUE): { OTHER_VAL:SCORE...}
 }`
 
 
 *compute empirical*
 
 * empirical[l] is a dictionary that maps a field value to the percentage of records with that value.
 
 `f ,counts = np.unique(self.x[:, l], return_counts=True)`
 
 