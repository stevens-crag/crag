
#include "ShortBraidForm.h"
#include "DehornoyForm.h"
#include "LinkedBraidStructure.h"
#include "ThRightNormalForm.h"
#include "braid_group.h"

LinkedBraidStructure shortenLBS(LinkedBraidStructure& lbs) {
  LinkedBraidStructure result = lbs;
  for (int i = 0; i < 4; ++i) {
    if (i % 2 == 0)
      lbs.removeRightHandles();
    else
      lbs.removeLeftHandles();
    if (result.size() > lbs.size())
      result = lbs;
  }
  return result;
}

Word dehornoy(int N, const Word& w) {
  LinkedBraidStructure df(N - 1, w);
  df.removeLeftHandles();
  return df.translateIntoWord();
}

Word garsideDehornoy(int N, const Word& w) {
  crag::braidgroup::BraidGroup B(N);
  ThRightNormalForm NF(B, w);
  const auto w1 = NF.getShortWord();
  return dehornoy(N, w1);
}

Word shortenBraid(int N, const Word& w) {
  LinkedBraidStructure df(N - 1, w);
  LinkedBraidStructure result = df;

  for (int i = 0; i < 4; ++i) {
    if (i % 2 == 0) {
      df.removeRightHandles();
    } else {
      df.removeLeftHandles();
    }

    if (df.size() < result.size()) {
      result = df;
    }
  }

  return result.translateIntoWord();
}

Word shortenBraid2(int n, const Word& w) {
  //  return shortenBraid(n, w);

  std::vector<int> result(w.begin(), w.end());

  for (size_t step = 16; step < result.size(); step *= 2) {
    std::vector<int> interm_result;
    interm_result.reserve(result.size());

    for (size_t i = 0; i < result.size(); i += step) {
      const auto begin = i;
      const auto end = std::min(i + step, result.size());

      Word seg(result.begin() + begin, result.begin() + end);

      const auto short_seg = shortenBraid(n, seg);

      if (seg.length() <= short_seg.length()) {
        interm_result.insert(interm_result.end(), result.begin() + begin, result.begin() + end);
      } else {
        interm_result.insert(interm_result.end(), short_seg.begin(), short_seg.end());
      }
    }

    if (interm_result.size() < result.size()) {
      result = interm_result;
    }
  }

  const auto result_word = Word(result);

  const auto short_word = shortenBraid(n, result_word);

  if (short_word.length() < result_word.length()) {
    return short_word;
  }

  return result_word;
}

Word shortBraidForm(int N, const Word& w) {
  crag::braidgroup::BraidGroup B(N);
  ThRightNormalForm NF(B, w);
  Word w1 = NF.getShortWord();
  //  cout << "int len = " << w1.length() << endl;
  return shortenBraid(N, w1);
}

vector<Word> shortBraidSbgpForm(int N, const vector<Word>& sbgp) {
  vector<Word> result(sbgp.size());
  for (int i = 0; i < sbgp.size(); ++i)
    result[i] = shortBraidForm(N, sbgp[i]);

  return result;
}
