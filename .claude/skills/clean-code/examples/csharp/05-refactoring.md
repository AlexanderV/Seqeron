# C# Refactoring Examples

> **Domain: Education Platform**
> Examples use Student, Course, Grade, Enrollment, Assignment, and Instructor entities.

## Code Smell: Long Method

**Before:**
```csharp
public void ProcessEnrollment(Enrollment enrollment)
{
    // 50+ lines of validation, calculations, notifications...
}
```

**After:**
```csharp
public void ProcessEnrollment(Enrollment enrollment)
{
    ValidateEnrollment(enrollment);
    CheckPrerequisites(enrollment);
    CalculateTuition(enrollment);
    NotifyStudent(enrollment);
}
```

---

## Code Smell: Duplicate Code

**Before:**
```csharp
public void EnrollStudent(Student student, Course course)
{
    if (string.IsNullOrEmpty(student.Email))
        throw new ValidationException("Email required");
    if (student.GPA < 2.0m)
        throw new ValidationException("GPA too low");
    // ...
}

public void TransferStudent(Student student, Course targetCourse)
{
    if (string.IsNullOrEmpty(student.Email))
        throw new ValidationException("Email required");
    if (student.GPA < 2.0m)
        throw new ValidationException("GPA too low");
    // ...
}
```

**After:**
```csharp
public void EnrollStudent(Student student, Course course)
{
    ValidateStudent(student);
    // ...
}

public void TransferStudent(Student student, Course targetCourse)
{
    ValidateStudent(student);
    // ...
}

private void ValidateStudent(Student student)
{
    if (string.IsNullOrEmpty(student.Email))
        throw new ValidationException("Email required");
    if (student.GPA < MinimumRequiredGPA)
        throw new ValidationException("GPA too low");
}
```

---

## Code Smell: Magic Numbers

**Before:**
```csharp
if (student.GPA > 3.7m)
    scholarshipAmount = tuition * 0.5m;

if (course.Credits > 4)
    maxEnrollment = 25;
```

**After:**
```csharp
private const decimal HonorsGPAThreshold = 3.7m;
private const decimal HonorsScholarshipRate = 0.5m;
private const int HighCreditThreshold = 4;
private const int ReducedClassSize = 25;

if (student.GPA > HonorsGPAThreshold)
    scholarshipAmount = tuition * HonorsScholarshipRate;

if (course.Credits > HighCreditThreshold)
    maxEnrollment = ReducedClassSize;
```

---

## Code Smell: Long Parameter List

**Before:**
```csharp
void CreateCourse(
    string code, 
    string title, 
    string description,
    int credits, 
    Guid instructorId, 
    string instructorName,
    string department, 
    string building, 
    string room,
    TimeSpan startTime, 
    TimeSpan endTime,
    List<DayOfWeek> meetingDays);
```

**After:**
```csharp
void CreateCourse(CourseCreationRequest request);

public record CourseCreationRequest
{
    public required CourseInfo Course { get; init; }
    public required InstructorInfo Instructor { get; init; }
    public required LocationInfo Location { get; init; }
    public required ScheduleInfo Schedule { get; init; }
}

public record CourseInfo(string Code, string Title, string Description, int Credits);
public record InstructorInfo(Guid Id, string Name, string Department);
public record LocationInfo(string Building, string Room);
public record ScheduleInfo(TimeSpan Start, TimeSpan End, List<DayOfWeek> Days);
```

---

## Code Smell: Switch Statement

**Before:**
```csharp
decimal CalculateGradePoints(Grade grade)
{
    switch (grade)
    {
        case Grade.A: return 4.0m;
        case Grade.B: return 3.0m;
        case Grade.C: return 2.0m;
        case Grade.D: return 1.0m;
        case Grade.F: return 0.0m;
        default: throw new ArgumentException();
    }
}
```

**After (Pattern Matching):**
```csharp
decimal CalculateGradePoints(Grade grade) => grade switch
{
    Grade.A or Grade.APlus => 4.0m,
    Grade.AMinus => 3.7m,
    Grade.BPlus => 3.3m,
    Grade.B => 3.0m,
    Grade.BMinus => 2.7m,
    Grade.CPlus => 2.3m,
    Grade.C => 2.0m,
    Grade.CMinus => 1.7m,
    Grade.D => 1.0m,
    Grade.F => 0.0m,
    _ => throw new UnknownGradeException(grade)
};
```

**After (Polymorphism for complex logic):**
```csharp
public abstract class AssessmentType
{
    public abstract decimal CalculateScore(Assignment assignment);
    public abstract bool AllowsLateSubmission { get; }
}

public class ExamAssessment : AssessmentType
{
    public override decimal CalculateScore(Assignment a) => a.PointsEarned / a.MaxPoints * 100;
    public override bool AllowsLateSubmission => false;
}

public class HomeworkAssessment : AssessmentType
{
    public override decimal CalculateScore(Assignment a) => 
        a.PointsEarned / a.MaxPoints * 100 * (a.IsLate ? 0.9m : 1.0m);
    public override bool AllowsLateSubmission => true;
}

public class ProjectAssessment : AssessmentType
{
    public override decimal CalculateScore(Assignment a) => 
        (a.PointsEarned + a.BonusPoints) / a.MaxPoints * 100;
    public override bool AllowsLateSubmission => true;
}
```

---

## Code Smell: Deep Nesting

**Before:**
```csharp
public void SubmitAssignment(Student student, Assignment assignment)
{
    if (student != null)
    {
        if (student.IsEnrolled)
        {
            if (assignment != null)
            {
                if (!assignment.IsPastDeadline())
                {
                    if (assignment.Content != null)
                    {
                        // Finally submit
                        _repository.Save(assignment);
                    }
                }
            }
        }
    }
}
```

**After (Guard Clauses):**
```csharp
public void SubmitAssignment(Student student, Assignment assignment)
{
    ArgumentNullException.ThrowIfNull(student);
    ArgumentNullException.ThrowIfNull(assignment);

    if (!student.IsEnrolled)
        throw new StudentNotEnrolledException(student.Id);

    if (assignment.IsPastDeadline())
        throw new DeadlinePassedException(assignment.Id, assignment.Deadline);

    if (assignment.Content is null)
        throw new EmptySubmissionException(assignment.Id);

    _repository.Save(assignment);
}
```

---

## Code Smell: Primitive Obsession

**Before:**
```csharp
public class Student
{
    public string Email { get; set; }           // Can be invalid
    public decimal GPA { get; set; }            // Can be negative
    public string StudentId { get; set; }       // Any string
}

// Validation scattered everywhere
if (!email.Contains("@")) throw new Exception();
if (gpa < 0 || gpa > 4) throw new Exception();
```

**After (Value Objects):**
```csharp
public record Email
{
    public string Value { get; }

    public Email(string value)
    {
        if (string.IsNullOrWhiteSpace(value) || !value.Contains('@'))
            throw new InvalidEmailException(value);
        Value = value.ToLowerInvariant();
    }

    public static implicit operator string(Email email) => email.Value;
}

public record GradePointAverage
{
    public decimal Value { get; }

    public GradePointAverage(decimal value)
    {
        if (value is < 0 or > 4.0m)
            throw new InvalidGPAException(value);
        Value = Math.Round(value, 2);
    }

    public bool IsHonors => Value >= 3.7m;
    public bool IsProbation => Value < 2.0m;
}

public record StudentId
{
    public string Value { get; }

    public StudentId(string value)
    {
        if (!Regex.IsMatch(value, @"^STU-\d{8}$"))
            throw new InvalidStudentIdException(value);
        Value = value;
    }
}

public class Student
{
    public StudentId Id { get; }
    public Email Email { get; }
    public GradePointAverage GPA { get; private set; }
}
```

---

## Code Smell: Feature Envy

**Before:**
```csharp
public class GradeCalculator
{
    public decimal CalculateFinalGrade(Student student)
    {
        // Accessing student's internal data too much
        var assignments = student.Assignments;
        var exams = student.Exams;
        var participation = student.ParticipationScore;

        var assignmentAvg = assignments.Average(a => a.Score);
        var examAvg = exams.Average(e => e.Score);

        return assignmentAvg * 0.3m + examAvg * 0.5m + participation * 0.2m;
    }
}
```

**After (Move method to data owner):**
```csharp
public class Student
{
    private readonly List<Assignment> _assignments = new();
    private readonly List<Exam> _exams = new();

    public decimal ParticipationScore { get; private set; }

    public decimal CalculateFinalGrade()
    {
        var assignmentAverage = CalculateAssignmentAverage();
        var examAverage = CalculateExamAverage();

        return assignmentAverage * AssignmentWeight
             + examAverage * ExamWeight
             + ParticipationScore * ParticipationWeight;
    }

    private decimal CalculateAssignmentAverage() =>
        _assignments.Count > 0 ? _assignments.Average(a => a.Score) : 0;

    private decimal CalculateExamAverage() =>
        _exams.Count > 0 ? _exams.Average(e => e.Score) : 0;

    private const decimal AssignmentWeight = 0.3m;
    private const decimal ExamWeight = 0.5m;
    private const decimal ParticipationWeight = 0.2m;
}
```

---

## Code Smell: Data Clumps

**Before:**
```csharp
public void ScheduleClass(
    string building, string room, int floor,      // Location clump
    int startHour, int startMinute,               // Time clump
    int endHour, int endMinute)
{
    // ...
}

public void ReserveRoom(
    string building, string room, int floor,      // Same clump again
    DateTime date)
{
    // ...
}
```

**After (Extract class):**
```csharp
public record Classroom(string Building, string Room, int Floor)
{
    public string FullLocation => $"{Building} {Room} (Floor {Floor})";
}

public record TimeSlot(TimeOnly Start, TimeOnly End)
{
    public TimeSpan Duration => End - Start;
    public bool OverlapsWith(TimeSlot other) =>
        Start < other.End && End > other.Start;
}

public void ScheduleClass(Classroom classroom, TimeSlot timeSlot)
{
    // Clean and focused
}

public void ReserveRoom(Classroom classroom, DateOnly date)
{
    // Reusing the same value object
}
```

---

## Refactoring Checklist

| Step | Action |
|------|--------|
| 1 | **Ensure tests exist** before refactoring |
| 2 | **Make small changes** one at a time |
| 3 | **Run tests** after each change |
| 4 | **Commit frequently** after successful refactorings |
| 5 | **Use IDE refactoring tools** (Rename, Extract Method, etc.) |

## Common Refactoring Patterns

| Smell | Refactoring |
|-------|-------------|
| Long Method | Extract Method |
| Duplicate Code | Extract Method, Pull Up Method |
| Magic Numbers | Extract Constant |
| Long Parameter List | Introduce Parameter Object |
| Switch Statement | Replace with Polymorphism |
| Deep Nesting | Guard Clauses, Extract Method |
| Primitive Obsession | Replace Primitive with Object |
| Feature Envy | Move Method |
| Data Clumps | Extract Class |

---

For architectural refactorings, see [Clean Architecture Patterns](../../../clean-architecture/PATTERNS.md).

For step-by-step refactoring journey, see [09-refactoring-journey.md](09-refactoring-journey.md).
