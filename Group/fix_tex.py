import re

file_path = r"D:\Apps\GItHubRebo\GeneralRelativityGraduationProject\test\Group.tex"
with open(file_path, "r", encoding="utf-8") as f:
    content = f.read()

# Fix sections
content = re.sub(r'\\section\*\{Section 1\.1:\s*(.*?)\}', r'\\chapter{Group Theory}\n\\section{\1}', content)
content = re.sub(r'\\section\*\{Section 2\.1:\s*(.*?)\}', r'\\chapter{Vector Spaces and Operators}\n\\section{\1}', content)
content = re.sub(r'\\section\*\{Section 3\.1:\s*(.*?)\}', r'\\chapter{Representation of Groups}\n\\section{\1}', content)
content = re.sub(r'\\section\*\{Section \d\.\d:\s*(.*?)\}', r'\\section{\1}', content)

# Function to properly format the mashed up Exam Answers
def format_exam_answer(text):
    # First, handle the title
    text = re.sub(r'^(Exam Answer:\s*.*?)(Statement of)', r'\\chapter*{\1}\n\\textbf{\2', text)
    # If the answer doesn't have "Statement of the Lemma" but has "Statement of the Theorem"
    text = text.replace('Statement of the Theorem:', '\\textbf{Statement of the Theorem:}\n')
    text = text.replace('Statement of the Lemma:', '\\textbf{Statement of the Lemma:}\n')
    text = text.replace('Proof:', '\n\\textbf{Proof:}\n')
    
    # Put Steps on new lines
    text = re.sub(r'(Step \d+(?: \([^)]+\))?:)', r'\n\1\n', text)
    
    # Handle $$ math blocks
    # We will split by $$ and alternate
    parts = text.split('$$')
    formatted_text = ""
    for i, part in enumerate(parts):
        if i % 2 == 1:
            formatted_text += '\n\\[\n' + part.strip() + '\n\\]\n'
        else:
            formatted_text += part
            
    # Fix some possible spacing issues
    formatted_text = re.sub(r'Case [A-Z]:', r'\n\\textbf{\g<0>}', formatted_text)
    formatted_text = formatted_text.replace('(Q.E.D.)', '\n\\textbf{(Q.E.D.)}\n')
    
    return formatted_text.strip()

# Find the lines that start with "Exam Answer:"
lines = content.split('\n')
new_lines = []
for line in lines:
    if line.startswith('Exam Answer:'):
        new_lines.append(format_exam_answer(line))
    else:
        new_lines.append(line)

content = '\n'.join(new_lines)

with open(file_path, "w", encoding="utf-8") as f:
    f.write(content)

print("Formatting applied successfully.")
