def generate_list(repeat_count):
    pattern = [30] + [1] * 30
    result = pattern * repeat_count
    return result

# Example usage:
repeat_count = 5
result_list = generate_list(repeat_count)
print(result_list)