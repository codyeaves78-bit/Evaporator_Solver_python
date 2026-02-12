Markdown
# Multiple Effect Evaporator Solver

This project simulates and optimizes multiple effect evaporator sets using **Harold Birkett's method**. It solves for juice distribution, brix profiles, pressure profiles, and heat transfer coefficients across multiple evaporator effects.
can handle 3 sets and 1 Pre

---

## 🏗 Installation

Follow these steps to set up the environment and run the simulation:

1. **Clone the Repository**
   ```bash
   git clone [https://github.com/your-username/evaporator-solver.git](https://github.com/your-username/evaporator-solver.git)
   cd evaporator-solver
Create a Virtual Environment (Recommended)

Bash
# Create the environment
python -m venv venv

# Activate on Windows:
venv\Scripts\activate
# Activate on macOS/Linux:
source venv/bin/activate
Install Dependencies
The project requires pandas, numpy, and scipy.

Bash
pip install -r requirements.txt
