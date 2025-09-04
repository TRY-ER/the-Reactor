import io
import logging
from functools import wraps
from rdkit import rdBase


def capture_rdkit_errors(func):
    """
    Decorator that captures RDKit error messages and ensures they're included
    in the return value when validation fails.
    
    This decorator works with functions that return tuple[bool, str] where:
    - First element is success/failure boolean
    - Second element is the message
    
    Usage:
        @capture_rdkit_errors
        def validate_something(self, input_str: str) -> tuple[bool, str]:
            # Your RDKit validation logic here
            return success, message
    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        # Setup logging to capture RDKit warnings/errors
        rdBase.LogToPythonLogger()
        
        # Create an in-memory stream and a handler to capture logs
        stream = io.StringIO()
        logger = logging.getLogger("rdkit")
        handler = logging.StreamHandler(stream)
        formatter = logging.Formatter("%(levelname)s: %(message)s")
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        logger.setLevel(logging.WARNING)  # Set the level to capture warnings
        
        try:
            result = func(*args, **kwargs)
            
            # Get any captured error messages
            captured_output = stream.getvalue().strip()
            
            # If function returned a tuple and validation failed
            if isinstance(result, tuple) and len(result) == 2:
                success, message = result
                if not success and captured_output:
                    # Enhance the error message with RDKit details
                    enhanced_message = f"{message}. RDKit Error: {captured_output}"
                    return False, enhanced_message
                else:
                    # Return original result
                    return result
            else:
                # Return original result if not the expected tuple format
                return result
                
        except Exception as e:
            captured_output = stream.getvalue().strip()
            if captured_output:
                return False, f"Exception: {str(e)}. RDKit Error: {captured_output}"
            else:
                return False, f"Exception: {str(e)}"
        finally:
            # Cleanup - remove the handler after we're done
            logger.removeHandler(handler)
    
    return wrapper
