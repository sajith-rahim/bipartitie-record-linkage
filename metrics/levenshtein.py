def levenshtein(self, str1, str2):
    if str1 is None:
        raise TypeError("str1 is None!")
    if str2 is None:
        raise TypeError("str2 is None!")
    if str1 == str2:
        return 0.0
    if len(str1) == 0:
        return len(str2)
    if len(str2) == 0:
        return len(str1)

    v0 = [0] * (len(str2) + 1)
    v1 = [0] * (len(str2) + 1)

    for i in range(len(v0)):
        v0[i] = i

    for i in range(len(str1)):
        v1[0] = i + 1
        for j in range(len(str2)):
            cost = 1
            if str1[i] == str2[j]:
                cost = 0
            v1[j + 1] = min(v1[j] + 1, v0[j + 1] + 1, v0[j] + cost)
        v0, v1 = v1, v0

    return v0[len(str2)]