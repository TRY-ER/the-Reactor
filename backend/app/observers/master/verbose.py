"""
This file implements the verbose mechanism for the observer to output relevant details to the user on run progress.
"""


def verbose_response(current_data, 
                     params_map, 
                     progress_map, 
                     execution_time=None):
    """
    This function prints the current status of the observation run. 
    It shows an intuitive progress update with a progress bar with colors. 

    You need to create a box. Within that box you need to have 3 progress section with bar.
    1. Model Type Progress
    2. Model Name Progress
    3. Reaction Type Progress

    At the bottom of each box the number of reactions found should be written.

    """
    reactions_found = len(current_data)

    model_name = current_data[0].get("model_name", "unknown")
    reaction_type = current_data[0].get("reaction_type", "unknown")

    model_type_index = progress_map.get("model_type_progress", 0)
    max_model_types = params_map.get("max_model_types", 0)

    model_name_index = progress_map.get("model_name_progress", 0)
    max_num_models = params_map.get("max_num_models", 0)

    reaction_type_index = progress_map.get("reaction_type_progress", 0)
    max_reaction_types = params_map.get("max_reaction_types", 67)

    def create_progress_bar(current, total, width=30):
        """Create a visual progress bar with colors"""
        if total == 0:
            percentage = 0
        else:
            percentage = min(current / total, 1.0)

        filled_width = int(percentage * width)
        empty_width = width - filled_width

        # Color codes
        green = '\033[92m'
        yellow = '\033[93m'
        red = '\033[91m'
        reset = '\033[0m'

        # Choose color based on progress
        if percentage >= 0.8:
            color = green
        elif percentage >= 0.4:
            color = yellow
        else:
            color = red

        bar = f"{color}{'█' * filled_width}{'░' * empty_width}{reset}"
        percent_text = f"{percentage * 100:.1f}%"
        return bar, percent_text

    # Create progress bars
    model_type_bar, model_type_percent = create_progress_bar(
        model_type_index, max_model_types)
    model_name_bar, model_name_percent = create_progress_bar(
        model_name_index, max_num_models)
    reaction_type_bar, reaction_type_percent = create_progress_bar(
        reaction_type_index, max_reaction_types)

    # Format execution time if provided
    time_display = ""
    if execution_time is not None and isinstance(execution_time, (int, float)):
        if execution_time < 60:
            time_display = f"Execution Time: {execution_time:.2f}s"
        else:
            minutes = int(execution_time // 60)
            seconds = execution_time % 60
            time_display = f"Execution Time: {minutes}m {seconds:.2f}s"
    else:
        execution_time = "N/A"

    print(f"""
┌───────────────────────────────────────────────────────────────┐
│                    EXPERIMENT PROGRESS                        │
├───────────────────────────────────────────────────────────────┤
  
  Model Name: {model_name} 
  Reaction Type: {reaction_type}  
  
  Model Type Progress:                                         
  {model_type_bar} {model_type_percent:>6}                     
  [{model_type_index:>3}/{max_model_types:<3}] Completed      
                                                               
  Model Name Progress:                                         
  {model_name_bar} {model_name_percent:>6}                     
  [{model_name_index:>3}/{max_num_models:<3}] Completed        
                                                               
  Reaction Type Progress:                                      
  {reaction_type_bar} {reaction_type_percent:>6}               
  [{reaction_type_index:>3}/{max_reaction_types:<3}] Completed
                                                               
├───────────────────────────────────────────────────────────────
│  Total Reactions Found: {reactions_found:>4}                  
│  {time_display:<59} 
└───────────────────────────────────────────────────────────────┘
    """)
