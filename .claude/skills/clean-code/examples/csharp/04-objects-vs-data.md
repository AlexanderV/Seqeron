# C# Objects vs Data Structures

Examples using **Healthcare** domain (Patient, Appointment, MedicalRecord).

For theory, see [Objects and Data Structures Principle](../../principles/05-objects-and-data-structures.md).

---

## Objects Hide Data, Expose Behavior

**❌ BAD - Anemic Domain Model:**
```csharp
public class Appointment
{
    public Guid Id { get; set; }
    public Guid PatientId { get; set; }
    public Guid DoctorId { get; set; }
    public DateTime ScheduledAt { get; set; }
    public AppointmentStatus Status { get; set; }
    public string? CancellationReason { get; set; }
}

// Business logic in service (wrong!)
public class AppointmentService
{
    public void Cancel(Appointment appointment, string reason)
    {
        if (appointment.Status == AppointmentStatus.Completed)
            throw new InvalidOperationException("Cannot cancel completed appointment");
        
        appointment.Status = AppointmentStatus.Cancelled;
        appointment.CancellationReason = reason;
    }
    
    public void Reschedule(Appointment appointment, DateTime newTime)
    {
        if (newTime < DateTime.UtcNow)
            throw new ArgumentException("Cannot schedule in the past");
        
        appointment.ScheduledAt = newTime;
        appointment.Status = AppointmentStatus.Rescheduled;
    }
}
```

**✅ GOOD - Rich Domain Model:**
```csharp
public class Appointment
{
    public AppointmentId Id { get; private set; }
    public PatientId PatientId { get; private set; }
    public DoctorId DoctorId { get; private set; }
    public DateTimeOffset ScheduledAt { get; private set; }
    public AppointmentStatus Status { get; private set; }
    public string? CancellationReason { get; private set; }
    
    private readonly List<DomainEvent> _domainEvents = [];
    
    private Appointment() { } // EF Core
    
    public static Appointment Schedule(
        PatientId patientId, 
        DoctorId doctorId, 
        DateTimeOffset scheduledAt)
    {
        if (scheduledAt < DateTimeOffset.UtcNow)
            throw new AppointmentSchedulingException("Cannot schedule in the past");
        
        var appointment = new Appointment
        {
            Id = AppointmentId.New(),
            PatientId = patientId,
            DoctorId = doctorId,
            ScheduledAt = scheduledAt,
            Status = AppointmentStatus.Scheduled
        };
        
        appointment._domainEvents.Add(new AppointmentScheduled(appointment.Id));
        return appointment;
    }
    
    public void Cancel(string reason)
    {
        if (Status == AppointmentStatus.Completed)
            throw new AppointmentStateException("Cannot cancel completed appointment");
        
        if (Status == AppointmentStatus.Cancelled)
            throw new AppointmentStateException("Appointment already cancelled");
        
        Status = AppointmentStatus.Cancelled;
        CancellationReason = reason;
        _domainEvents.Add(new AppointmentCancelled(Id, reason));
    }
    
    public void Reschedule(DateTimeOffset newTime)
    {
        if (Status != AppointmentStatus.Scheduled)
            throw new AppointmentStateException("Can only reschedule scheduled appointments");
        
        if (newTime < DateTimeOffset.UtcNow)
            throw new AppointmentSchedulingException("Cannot reschedule to past");
        
        var previousTime = ScheduledAt;
        ScheduledAt = newTime;
        _domainEvents.Add(new AppointmentRescheduled(Id, previousTime, newTime));
    }
    
    public void MarkAsCompleted()
    {
        if (Status != AppointmentStatus.Scheduled)
            throw new AppointmentStateException("Can only complete scheduled appointments");
        
        Status = AppointmentStatus.Completed;
        _domainEvents.Add(new AppointmentCompleted(Id));
    }
}
```

---

## Data Structures Expose Data, No Behavior

**✅ GOOD - DTOs for API:**
```csharp
// Request DTO - data structure for API input
public record ScheduleAppointmentRequest(
    Guid PatientId,
    Guid DoctorId,
    DateTime ScheduledAt);

// Response DTO - data structure for API output
public record AppointmentResponse(
    Guid Id,
    Guid PatientId,
    string PatientName,
    Guid DoctorId,
    string DoctorName,
    DateTime ScheduledAt,
    string Status);

// Mapping in handler
public class ScheduleAppointmentHandler(
    IPatientRepository patients,
    IDoctorRepository doctors,
    IAppointmentRepository appointments,
    IUnitOfWork unitOfWork)
{
    public async Task<AppointmentResponse> Handle(ScheduleAppointmentRequest request)
    {
        var patient = await patients.GetByIdAsync(new PatientId(request.PatientId))
            ?? throw new EntityNotFoundException<PatientId>("Patient", new PatientId(request.PatientId));
        
        var doctor = await doctors.GetByIdAsync(new DoctorId(request.DoctorId))
            ?? throw new EntityNotFoundException<DoctorId>("Doctor", new DoctorId(request.DoctorId));
        
        var appointment = Appointment.Schedule(
            patient.Id,
            doctor.Id,
            request.ScheduledAt);
        
        await appointments.AddAsync(appointment);
        await unitOfWork.SaveChangesAsync();
        
        return new AppointmentResponse(
            appointment.Id.Value,
            patient.Id.Value,
            patient.FullName,
            doctor.Id.Value,
            doctor.FullName,
            appointment.ScheduledAt.DateTime,
            appointment.Status.ToString());
    }
}
```

---

## Patient as Rich Domain Entity

```csharp
public class Patient
{
    public PatientId Id { get; private set; }
    public string FirstName { get; private set; }
    public string LastName { get; private set; }
    public DateOnly DateOfBirth { get; private set; }
    public Email Email { get; private set; }
    public PhoneNumber Phone { get; private set; }
    
    private readonly List<MedicalRecord> _medicalRecords = [];
    public IReadOnlyList<MedicalRecord> MedicalRecords => _medicalRecords.AsReadOnly();
    
    public string FullName => $"{FirstName} {LastName}";
    public int Age => CalculateAge(DateOfBirth);
    
    public void AddMedicalRecord(MedicalRecord record)
    {
        ArgumentNullException.ThrowIfNull(record);
        
        if (record.PatientId != Id)
            throw new InvalidOperationException("Record belongs to different patient");
        
        _medicalRecords.Add(record);
    }
    
    public void UpdateContactInfo(Email email, PhoneNumber phone)
    {
        Email = email ?? throw new ArgumentNullException(nameof(email));
        Phone = phone ?? throw new ArgumentNullException(nameof(phone));
    }
    
    private static int CalculateAge(DateOnly birthDate)
    {
        var today = DateOnly.FromDateTime(DateTime.Today);
        var age = today.Year - birthDate.Year;
        if (birthDate > today.AddYears(-age)) age--;
        return age;
    }
}
```

---

## Value Objects for Type Safety

```csharp
public readonly record struct PatientId(Guid Value)
{
    public static PatientId New() => new(Guid.NewGuid());
    public override string ToString() => Value.ToString()[..8];
}

public readonly record struct DoctorId(Guid Value)
{
    public static DoctorId New() => new(Guid.NewGuid());
}

public readonly record struct AppointmentId(Guid Value)
{
    public static AppointmentId New() => new(Guid.NewGuid());
}

public record Email
{
    public string Value { get; }
    
    public Email(string value)
    {
        if (string.IsNullOrWhiteSpace(value))
            throw new ArgumentException("Email cannot be empty");
        
        if (!value.Contains('@'))
            throw new ArgumentException("Invalid email format");
        
        Value = value.ToLowerInvariant();
    }
    
    public override string ToString() => Value;
}

public record PhoneNumber
{
    public string Value { get; }
    
    public PhoneNumber(string value)
    {
        var digits = new string(value.Where(char.IsDigit).ToArray());
        
        if (digits.Length < 10)
            throw new ArgumentException("Phone number must have at least 10 digits");
        
        Value = digits;
    }
    
    public string Formatted => $"({Value[..3]}) {Value[3..6]}-{Value[6..]}";
}
```

---

## Law of Demeter

**❌ BAD - Train Wreck:**
```csharp
string diagnosis = appointment.GetPatient().GetMedicalRecords().Last().GetDiagnosis().GetName();
```

**✅ GOOD - Tell, Don't Ask:**
```csharp
// Ask for what you need
string diagnosis = appointment.GetLastDiagnosis();

// In Appointment class
public string? GetLastDiagnosis()
{
    // Appointment delegates to its collaborators
    return _patient.GetLatestDiagnosis();
}

// In Patient class
public string? GetLatestDiagnosis()
{
    return _medicalRecords
        .OrderByDescending(r => r.CreatedAt)
        .FirstOrDefault()
        ?.Diagnosis;
}
```

---

## Summary

| Concept | Objects | Data Structures |
|---------|---------|-----------------|
| **Example** | `Appointment`, `Patient` | `AppointmentResponse`, `ScheduleAppointmentRequest` |
| **Data** | Private/protected | Public |
| **Behavior** | Public methods | None (or mapping only) |
| **Invariants** | Enforced internally | None |
| **Use case** | Domain logic | Data transfer, serialization |

## Related

- [Objects and Data Structures](../../principles/05-objects-and-data-structures.md) — Theory
- [Complete Example](08-complete-example.md) — Banking domain integration
- [Clean Architecture - DDD](../../../clean-architecture/principles/03-domain-driven-design.md) — Architectural patterns
