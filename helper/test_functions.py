from functions import rescale_threshold_if_out_of_range, is_zero_speciation_events_failure


def test_threshold_within_range_passes_through_unchanged():
    assert rescale_threshold_if_out_of_range(50, (0, 100)) == 50


def test_threshold_exactly_at_max_is_not_rescaled():
    assert rescale_threshold_if_out_of_range(100, (0, 100)) == 100


def test_threshold_over_max_on_fractional_scale_tree_is_rescaled():
    # e.g. FastTree's own local-support scale is natively 0-1, not 0-100.
    assert rescale_threshold_if_out_of_range(50, (0, 1)) == 0.5


def test_zero_threshold_is_never_rescaled():
    assert rescale_threshold_if_out_of_range(0, (0, 1)) == 0


def test_detects_the_known_zero_speciation_events_possvm_crash():
    log = "2026-09-17 [INFO] There are no speciation events in this tree.\nTraceback...\nKeyError: 'in_gene'"
    assert is_zero_speciation_events_failure(log) is True


def test_does_not_flag_an_unrelated_failure():
    log = "Traceback...\nFileNotFoundError: [Errno 2] No such file or directory: 'x.fasta'"
    assert is_zero_speciation_events_failure(log) is False
