void gsort(int *start, int *end)
{
	const int n = end - start;

	int i = 1;
	while (i < n) {
		if (start[i] < start[i-1]) {
			int tmp = start[i];
			start[i] = start[i-1];
			start[i-1] = tmp;
			if (i > 1) {
				i--;
			}
		} else {
			i++;
		}
	}
}
