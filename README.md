# NE_5001 Project

To run this project, please follow the instructions below.
## Prerequisites
-  Geant4 installation
-  CMake
-  A C++ compiler (e.g., g++, clang++)
## Building the Project
1. Clone the repository
2. Create a build directory:
   ```bash
   mkdir build
   cd build
   ```
3. Run CMake to configure the project:
    ```bash
    cmake ..
    ```
4. Build the project:
    ```bash
    make
    ```

## How to use git:
1. To clone the repository:
   ```bash
   git clone <repository_url>
   ```
2. To check the status of your repository:
   ```bash
   git status
   ```
3. To add changes to staging:
   ```bash
   git add <file_name>
   ```
   2. To add all changes to staging:
        ```bash
        git add .
        ```

4. To commit changes:
   ```bash
   git commit -m "Your commit message"
   ```
5. To push changes to the remote repository:
```bash
   git push
   ```
6. To pull the latest changes from the remote repository:
   ```bash
   git pull
   ```

To update the code in your local repository with the latest changes from the remote repository, use the `git pull` command. This command fetches the latest changes and merges them into your local branch.

To push local changes into the remote repository, follow steps 3-5.

In order for git to work you should be placed in the root directory of your project where the .git folder is located. In this case the "NE_5001_Project" folder.

Good practices when using git:
- Commit often with clear messages.
- Pull before pushing to avoid conflicts.
