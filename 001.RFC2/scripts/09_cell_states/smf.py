"""Use an already prepared scMalignantFinder classifier; setup and file I/O are caller responsibilities."""

def predict_malignancy(prepared_classifier):
    """Return the predicted AnnData object without loading or writing files here."""
    return prepared_classifier.predict()
