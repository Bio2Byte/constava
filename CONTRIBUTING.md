# Contributing to Constava

Contributions of all kinds are welcome: bug reports, documentation improvements, feature requests, and code changes.

Please read the guidelines below to make the contribution process smooth for everyone.

---

## How to Contribute

### Reporting Bugs

If you find a bug, please open an issue and include:
- A clear and descriptive title
- Steps to reproduce the issue
- Expected and actual behavior
- Your environment (OS, Python version, Constava version)

If possible, include a minimal example or traceback.

### Suggesting Enhancements

Feature requests and improvements are welcome.  
When opening an issue, try to explain:
- The problem you are trying to solve
- Why the current behavior is insufficient
- Any alternative solutions you have considered

### Development Setup

1. Fork the repository and clone your fork.
2. Create a virtual environment.
3. Install dependencies:
   ```bash
   pip install -r requirements.txt
   pip install -e .
   ```

4.	Run the test suite to ensure everything works:
    ```bash
    constava test
    ```

### Making Changes

1. Create a new branch for your changes.
1. Keep commits focused and atomic.
1. Write or update tests when changing behavior.
1. Update documentation if your change affects public APIs or CLI usage.


### Code Style

1. Follow standard Python conventions (PEP 8).
1. Keep code readable and well-documented.
1. Avoid unnecessary complexity; clarity beats cleverness.

### Testing

All changes should pass the test suite before submission:

```bash
constava test
```

Pull requests that break tests will not be merged.


### Submitting a Pull Request

1. Open a pull request against the default branch.
1. Clearly describe what the PR does and why it is needed.
1. Link related issues where applicable.
1. Be prepared to revise your PR based on review feedback.


## Code of Conduct

This project follows a Code of Conduct.
By participating, you are expected to uphold it in all interactions.
