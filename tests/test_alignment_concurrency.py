# Regression test for concurrent alignments sharing one set of temp files.

import threading

from services.alignment import align_with_mafft

# Two inputs that are easy to tell apart: if a call returns the other one's alleles, or
# the single blank record MAFFT's empty output parses into, the two runs collided.
FIRST = {"first_a": "ACGTACGTACGTAAAA", "first_b": "ACGTACGTACGTAAAC"}
SECOND = {"second_a": "TTTTGGGGCCCCTTTT", "second_b": "TTTTGGGGCCCCTTTA",
          "second_c": "TTTTGGGGCCCCTTAA"}

def test_concurrent_alignments_do_not_read_each_others_files():
    """Two alignments at once must each get their own answer.

    align_with_mafft used to write 'tmp/unaligned_fasta_tmp.fasta' and read
    'tmp/aligned_fasta_tmp.fasta' - fixed names under the working directory, which every
    gunicorn worker shares. Run concurrently, one call overwrote the other's input and read
    its output: the loser returned {'': ''} or the other request's alignment, as a 200.

    Exceptions are collected rather than left to propagate, because an exception in a thread
    does not fail the test that started it - it is printed and the thread ends. The old code
    raises here as well as answering wrongly (two callers racing to create the directory),
    so a version of this test that only compared return values passed against it.
    """
    problems = []

    def align_repeatedly(sequences):
        expected = sorted(sequences)
        for _ in range(15):
            try:
                result = align_with_mafft(sequences)
            except Exception as error:
                problems.append(f"{expected[0]}: raised {type(error).__name__}")
                continue
            if sorted(result) != expected:
                problems.append(f"{expected[0]}: expected {expected}, got {sorted(result)}")

    threads = [threading.Thread(target=align_repeatedly, args=(sequences,))
               for sequences in (FIRST, SECOND)]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()

    assert problems == []
