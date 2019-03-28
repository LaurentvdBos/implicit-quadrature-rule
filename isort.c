void isort(int *start, int *end)
{
	int n = end - start;

	for (int i = 1; i < n; i++) {
		for (int j = i; j > 0 && start[j] < start[j-1]; j--) {
			int tmp = start[j];
			start[j] = start[j-1];
			start[j-1] = tmp;
		}
	}
}

