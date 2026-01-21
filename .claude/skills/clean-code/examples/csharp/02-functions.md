# C# Functions Examples

> **Domain: IoT Sensor Platform**
> Examples use Device, Sensor, Reading, Alert, Telemetry, and Threshold entities.

## Table of Contents

1. [Small Functions](#small-functions)
2. [Do One Thing](#do-one-thing)
3. [One Level of Abstraction](#one-level-of-abstraction)
4. [Function Arguments](#function-arguments)
5. [Command-Query Separation](#command-query-separation)
6. [Avoid Side Effects](#avoid-side-effects)
7. [Prefer Exceptions to Error Codes](#prefer-exceptions-to-error-codes)
8. [Don't Repeat Yourself](#dont-repeat-yourself)

---

## Small Functions

Functions should be small - really small. 5-20 lines is ideal.

### Example 1: Sensor Reading Processing

**❌ BAD - Long Function:**
```csharp
public async Task ProcessSensorReading(SensorReading reading)
{
    // Validate
    if (reading == null)
        throw new ArgumentNullException(nameof(reading));

    if (reading.SensorId == Guid.Empty)
        throw new InvalidOperationException("Reading must have a sensor ID");

    if (reading.Value < -1000 || reading.Value > 1000)
        throw new InvalidOperationException("Reading value out of range");

    // Get sensor configuration
    var sensor = await _sensorRepository.GetByIdAsync(reading.SensorId);
    if (sensor == null)
        throw new InvalidOperationException($"Sensor {reading.SensorId} not found");

    if (!sensor.IsActive)
        throw new InvalidOperationException($"Sensor {reading.SensorId} is inactive");

    // Check thresholds
    var threshold = await _thresholdRepository.GetBySensorTypeAsync(sensor.Type);
    bool isAnomalous = false;

    if (reading.Value > threshold.MaxValue)
    {
        isAnomalous = true;
        var alert = new Alert
        {
            SensorId = sensor.Id,
            Type = AlertType.HighValue,
            Message = $"Value {reading.Value} exceeds maximum {threshold.MaxValue}",
            Timestamp = DateTime.UtcNow
        };
        await _alertRepository.AddAsync(alert);
        await _notificationService.SendAsync(sensor.OwnerId, alert.Message);
    }

    if (reading.Value < threshold.MinValue)
    {
        isAnomalous = true;
        var alert = new Alert
        {
            SensorId = sensor.Id,
            Type = AlertType.LowValue,
            Message = $"Value {reading.Value} below minimum {threshold.MinValue}",
            Timestamp = DateTime.UtcNow
        };
        await _alertRepository.AddAsync(alert);
        await _notificationService.SendAsync(sensor.OwnerId, alert.Message);
    }

    // Store reading
    reading.IsAnomalous = isAnomalous;
    await _readingRepository.AddAsync(reading);

    // Update sensor stats
    sensor.LastReadingAt = DateTime.UtcNow;
    sensor.TotalReadings++;
    await _sensorRepository.UpdateAsync(sensor);
}
```

**✅ GOOD - Small Functions:**
```csharp
public async Task ProcessSensorReading(SensorReading reading)
{
    ValidateReading(reading);

    var sensor = await GetActiveSensorAsync(reading.SensorId);
    var threshold = await GetThresholdAsync(sensor.Type);

    var isAnomalous = await CheckAndAlertThresholdViolations(reading, sensor, threshold);

    await StoreReadingAsync(reading, isAnomalous);
    await UpdateSensorStatsAsync(sensor);
}

private void ValidateReading(SensorReading reading)
{
    ArgumentNullException.ThrowIfNull(reading);

    if (reading.SensorId == Guid.Empty)
        throw new InvalidReadingException("Reading must have a sensor ID");

    if (reading.Value is < MinValidValue or > MaxValidValue)
        throw new InvalidReadingException($"Reading value {reading.Value} out of valid range");
}

private async Task<Sensor> GetActiveSensorAsync(SensorId sensorId)
{
    var sensor = await _sensorRepository.GetByIdAsync(sensorId)
        ?? throw new SensorNotFoundException(sensorId);

    if (!sensor.IsActive)
        throw new InactiveSensorException(sensorId);

    return sensor;
}

private async Task<Threshold> GetThresholdAsync(SensorType sensorType)
{
    return await _thresholdRepository.GetBySensorTypeAsync(sensorType);
}

private async Task<bool> CheckAndAlertThresholdViolations(
    SensorReading reading, 
    Sensor sensor, 
    Threshold threshold)
{
    bool isAnomalous = false;

    if (reading.Value > threshold.MaxValue)
    {
        isAnomalous = true;
        await CreateAndSendAlertAsync(sensor, AlertType.HighValue, 
            $"Value {reading.Value} exceeds maximum {threshold.MaxValue}");
    }

    if (reading.Value < threshold.MinValue)
    {
        isAnomalous = true;
        await CreateAndSendAlertAsync(sensor, AlertType.LowValue, 
            $"Value {reading.Value} below minimum {threshold.MinValue}");
    }

    return isAnomalous;
}

private async Task CreateAndSendAlertAsync(Sensor sensor, AlertType type, string message)
{
    var alert = new Alert(sensor.Id, type, message);
    await _alertRepository.AddAsync(alert);
    await _notificationService.SendAsync(sensor.OwnerId, message);
}

private async Task StoreReadingAsync(SensorReading reading, bool isAnomalous)
{
    reading.IsAnomalous = isAnomalous;
    await _readingRepository.AddAsync(reading);
}

private async Task UpdateSensorStatsAsync(Sensor sensor)
{
    sensor.RecordReading();
    await _sensorRepository.UpdateAsync(sensor);
}
```

---

## Do One Thing

Functions should do one thing, do it well, and do it only.

### Example 2: Device Registration

**❌ BAD - Does Multiple Things:**
```csharp
public async Task<bool> RegisterDevice(string serialNumber, string firmwareVersion, DeviceType type)
{
    // Check if device exists
    var existingDevice = await _deviceRepository.GetBySerialAsync(serialNumber);
    if (existingDevice != null)
        return false;

    // Validate serial
    if (serialNumber.Length != 12)
        return false;

    // Validate firmware
    if (!Version.TryParse(firmwareVersion, out var version))
        return false;

    // Check firmware compatibility
    var minVersion = await _firmwareService.GetMinVersionAsync(type);
    if (version < minVersion)
        return false;

    // Create device
    var device = new Device
    {
        SerialNumber = serialNumber,
        FirmwareVersion = firmwareVersion,
        Type = type,
        RegisteredAt = DateTime.UtcNow,
        Status = DeviceStatus.Pending
    };

    // Save to database
    await _deviceRepository.AddAsync(device);

    // Provision cloud resources
    await _cloudService.ProvisionAsync(device.Id);

    // Generate credentials
    var credentials = await _authService.GenerateDeviceCredentialsAsync(device.Id);
    await _credentialRepository.AddAsync(credentials);

    // Log
    _logger.LogInformation($"Device registered: {serialNumber}");

    return true;
}
```

**✅ GOOD - Each Function Does One Thing:**
```csharp
public async Task RegisterDevice(DeviceRegistrationRequest request)
{
    await EnsureSerialNumberIsUnique(request.SerialNumber);
    ValidateDeviceInput(request);
    await ValidateFirmwareCompatibility(request);

    var device = CreateDevice(request);
    await SaveDeviceAsync(device);
    await ProvisionCloudResourcesAsync(device);
    await GenerateDeviceCredentialsAsync(device);

    LogDeviceRegistration(device);
}

private async Task EnsureSerialNumberIsUnique(string serialNumber)
{
    var existing = await _deviceRepository.GetBySerialAsync(serialNumber);
    if (existing != null)
        throw new DuplicateSerialNumberException(serialNumber);
}

private void ValidateDeviceInput(DeviceRegistrationRequest request)
{
    if (request.SerialNumber.Length != SerialNumberLength)
        throw new InvalidSerialNumberException(request.SerialNumber);

    if (!Version.TryParse(request.FirmwareVersion, out _))
        throw new InvalidFirmwareVersionException(request.FirmwareVersion);
}

private async Task ValidateFirmwareCompatibility(DeviceRegistrationRequest request)
{
    var minVersion = await _firmwareService.GetMinVersionAsync(request.Type);
    var deviceVersion = Version.Parse(request.FirmwareVersion);

    if (deviceVersion < minVersion)
        throw new IncompatibleFirmwareException(request.FirmwareVersion, minVersion);
}

private Device CreateDevice(DeviceRegistrationRequest request)
{
    return new Device
    {
        SerialNumber = request.SerialNumber,
        FirmwareVersion = request.FirmwareVersion,
        Type = request.Type,
        RegisteredAt = DateTime.UtcNow,
        Status = DeviceStatus.Pending
    };
}
```

---

## One Level of Abstraction

Don't mix high-level and low-level operations in the same function.

### Example 3: Telemetry Report Generation

**❌ BAD - Mixed Abstraction Levels:**
```csharp
public string GenerateDeviceTelemetryReport(DeviceId deviceId, DateRange range)
{
    // High-level concept
    var readings = _readingRepository.GetByDeviceAndRange(deviceId, range);

    // Low-level JSON building
    var json = "{\"device_id\":\"" + deviceId + "\",";
    json += "\"readings\":[";

    // Medium-level calculation
    var avgTemperature = readings.Where(r => r.Type == SensorType.Temperature).Average(r => r.Value);

    // Low-level string manipulation
    json += "{\"avg_temperature\":" + avgTemperature.ToString("F2") + ",";
    json += "\"readings_count\":" + readings.Count + ",";
    json += "\"data\":[";

    foreach (var reading in readings)
    {
        json += "{\"timestamp\":\"" + reading.Timestamp.ToString("O") + "\",";
        json += "\"value\":" + reading.Value + ",";
        json += "\"sensor\":\"" + reading.SensorType + "\"},";
    }

    json = json.TrimEnd(',') + "]}]}";
    return json;
}
```

**✅ GOOD - Consistent Abstraction Level:**
```csharp
// High level - orchestration only
public string GenerateDeviceTelemetryReport(DeviceId deviceId, DateRange range)
{
    var telemetryData = CollectTelemetryData(deviceId, range);
    var report = BuildTelemetryReport(telemetryData);
    return SerializeReport(report);
}

// Medium level - data collection
private TelemetryData CollectTelemetryData(DeviceId deviceId, DateRange range)
{
    var readings = _readingRepository.GetByDeviceAndRange(deviceId, range);

    return new TelemetryData
    {
        DeviceId = deviceId,
        Range = range,
        Readings = readings,
        Statistics = CalculateStatistics(readings)
    };
}

private TelemetryStatistics CalculateStatistics(IReadOnlyList<SensorReading> readings)
{
    return new TelemetryStatistics
    {
        AverageTemperature = CalculateAverageByType(readings, SensorType.Temperature),
        AverageHumidity = CalculateAverageByType(readings, SensorType.Humidity),
        TotalReadings = readings.Count
    };
}

private double? CalculateAverageByType(IReadOnlyList<SensorReading> readings, SensorType type)
{
    var filtered = readings.Where(r => r.Type == type).ToList();
    return filtered.Count > 0 ? filtered.Average(r => r.Value) : null;
}

// Medium level - report building
private TelemetryReport BuildTelemetryReport(TelemetryData data)
{
    return new TelemetryReport
    {
        DeviceId = data.DeviceId,
        GeneratedAt = DateTime.UtcNow,
        Statistics = data.Statistics,
        Readings = data.Readings.Select(MapToReportReading).ToList()
    };
}

private ReportReading MapToReportReading(SensorReading reading)
{
    return new ReportReading(reading.Timestamp, reading.Value, reading.Type);
}

// Low level - serialization
private string SerializeReport(TelemetryReport report)
{
    return JsonSerializer.Serialize(report, _jsonOptions);
}
```

---

## Function Arguments

Ideal: 0-1 arguments. Acceptable: 2. Avoid: 3+.

### Example 4: Parameter Objects

**❌ BAD - Too Many Parameters:**
```csharp
public void ConfigureSensor(
    Guid sensorId,
    string name,
    SensorType type,
    double minThreshold,
    double maxThreshold,
    int samplingRateMs,
    bool enableAlerts,
    string alertEmail,
    int alertCooldownMinutes,
    bool enableDataLogging,
    int dataRetentionDays)
{
    // ...
}
```

**✅ GOOD - Parameter Object:**
```csharp
public void ConfigureSensor(SensorConfiguration config)
{
    ValidateConfiguration(config);
    ApplyConfiguration(config);
}

public record SensorConfiguration
{
    public required SensorId SensorId { get; init; }
    public required string Name { get; init; }
    public required SensorType Type { get; init; }
    public required ThresholdSettings Thresholds { get; init; }
    public required SamplingSettings Sampling { get; init; }
    public required AlertSettings Alerts { get; init; }
    public required DataLoggingSettings DataLogging { get; init; }
}

public record ThresholdSettings(double MinValue, double MaxValue);
public record SamplingSettings(int RateMilliseconds);
public record AlertSettings(bool Enabled, string Email, int CooldownMinutes);
public record DataLoggingSettings(bool Enabled, int RetentionDays);

// Clean usage with object initializer
var config = new SensorConfiguration
{
    SensorId = sensorId,
    Name = "Temperature Sensor A1",
    Type = SensorType.Temperature,
    Thresholds = new ThresholdSettings(-40, 85),
    Sampling = new SamplingSettings(1000),
    Alerts = new AlertSettings(true, "ops@company.com", 15),
    DataLogging = new DataLoggingSettings(true, 90)
};

ConfigureSensor(config);
```

### Example 5: Avoid Flag Arguments

**❌ BAD - Boolean Flag:**
```csharp
public async Task SendReading(SensorReading reading, bool immediate)
{
    if (immediate)
    {
        await _messageQueue.SendImmediateAsync(reading);
    }
    else
    {
        _messageQueue.Enqueue(reading);
    }
}

// Confusing usage
await SendReading(reading, true); // What does true mean?
```

**✅ GOOD - Separate Methods:**
```csharp
public async Task SendReadingImmediately(SensorReading reading)
{
    await _messageQueue.SendImmediateAsync(reading);
}

public void EnqueueReading(SensorReading reading)
{
    _messageQueue.Enqueue(reading);
}

// Clear usage
await SendReadingImmediately(criticalReading);
EnqueueReading(normalReading);
```

---

## Command-Query Separation

Functions should either do something or answer something, not both.

### Example 6: CQS Applied to Device Status

**❌ BAD - Mixed Command and Query:**
```csharp
// Returns status AND updates last check time - confusing!
public DeviceStatus CheckDeviceStatus(DeviceId deviceId)
{
    var device = _repository.GetById(deviceId);
    device.LastCheckedAt = DateTime.UtcNow;  // Side effect!
    _repository.Update(device);
    return device.Status;
}
```

**✅ GOOD - Separated:**
```csharp
// Query - returns information, no side effects
public DeviceStatus GetDeviceStatus(DeviceId deviceId)
{
    var device = _repository.GetById(deviceId);
    return device.Status;
}

// Query - check connectivity
public async Task<bool> IsDeviceOnline(DeviceId deviceId)
{
    var lastPing = await _pingService.GetLastPingAsync(deviceId);
    return lastPing > DateTime.UtcNow.AddMinutes(-5);
}

// Command - explicitly updates state
public async Task UpdateLastHealthCheck(DeviceId deviceId)
{
    var device = await _repository.GetByIdAsync(deviceId);
    device.LastCheckedAt = DateTime.UtcNow;
    await _repository.UpdateAsync(device);
}

// Clear usage
var status = GetDeviceStatus(deviceId);
if (status == DeviceStatus.Active && await IsDeviceOnline(deviceId))
{
    await UpdateLastHealthCheck(deviceId);
}
```

---

## Avoid Side Effects

Functions shouldn't have hidden side effects.

### Example 7: Hidden Side Effects in Sensor Validation

**❌ BAD - Hidden Side Effect:**
```csharp
public bool ValidateSensorConnection(SensorId sensorId)
{
    var sensor = _repository.GetById(sensorId);

    if (sensor.TestConnection())
    {
        // HIDDEN SIDE EFFECT! Name says "validate" but also calibrates!
        sensor.Calibrate();
        sensor.Status = SensorStatus.Ready;
        _repository.Update(sensor);
        return true;
    }

    return false;
}
```

**✅ GOOD - Explicit:**
```csharp
public bool ValidateSensorConnection(SensorId sensorId)
{
    var sensor = _repository.GetById(sensorId);
    return sensor.TestConnection();
}

public async Task InitializeSensorAsync(SensorId sensorId)
{
    if (!ValidateSensorConnection(sensorId))
        throw new SensorConnectionException(sensorId);

    var sensor = await _repository.GetByIdAsync(sensorId);
    await CalibrateSensorAsync(sensor);
    await ActivateSensorAsync(sensor);
}

private async Task CalibrateSensorAsync(Sensor sensor)
{
    await sensor.CalibrateAsync();
    _logger.LogInformation("Sensor {Id} calibrated", sensor.Id);
}

private async Task ActivateSensorAsync(Sensor sensor)
{
    sensor.Status = SensorStatus.Ready;
    await _repository.UpdateAsync(sensor);
}
```

---

## Prefer Exceptions to Error Codes

**❌ BAD - Error Codes:**
```csharp
public const int SUCCESS = 0;
public const int ERROR_SENSOR_NOT_FOUND = 1;
public const int ERROR_SENSOR_OFFLINE = 2;
public const int ERROR_THRESHOLD_EXCEEDED = 3;

public int TakeReading(SensorId sensorId)
{
    var sensor = _repository.GetById(sensorId);
    if (sensor == null) return ERROR_SENSOR_NOT_FOUND;
    if (!sensor.IsOnline) return ERROR_SENSOR_OFFLINE;

    var value = sensor.Read();
    if (value > sensor.MaxThreshold) return ERROR_THRESHOLD_EXCEEDED;

    return SUCCESS;
}
```

**✅ GOOD - Exceptions with Domain Types:**
```csharp
public SensorReading TakeReading(SensorId sensorId)
{
    var sensor = GetOnlineSensorOrThrow(sensorId);
    var reading = sensor.TakeReading();
    ValidateReading(reading, sensor.Threshold);
    return reading;
}

private Sensor GetOnlineSensorOrThrow(SensorId sensorId)
{
    var sensor = _repository.GetById(sensorId)
        ?? throw new SensorNotFoundException(sensorId);

    if (!sensor.IsOnline)
        throw new SensorOfflineException(sensorId);

    return sensor;
}

private void ValidateReading(SensorReading reading, Threshold threshold)
{
    if (reading.Value > threshold.MaxValue)
        throw new ThresholdExceededException(reading, threshold);
}
```

---

## Don't Repeat Yourself

### Example 8: Repeated Alert Logic

**❌ BAD - Duplicated Alert Creation:**
```csharp
public async Task CheckTemperatureSensor(SensorId sensorId)
{
    var reading = await GetLatestReadingAsync(sensorId);
    if (reading.Value > 80)
    {
        var alert = new Alert
        {
            SensorId = sensorId,
            Type = AlertType.Critical,
            Message = "Temperature exceeds safe limit",
            Timestamp = DateTime.UtcNow,
            Priority = Priority.High
        };
        await _alertRepository.AddAsync(alert);
        await _notificationService.NotifyAsync(alert);
    }
}

public async Task CheckHumiditySensor(SensorId sensorId)
{
    var reading = await GetLatestReadingAsync(sensorId);
    if (reading.Value > 90)
    {
        var alert = new Alert
        {
            SensorId = sensorId,
            Type = AlertType.Critical,
            Message = "Humidity exceeds safe limit",
            Timestamp = DateTime.UtcNow,
            Priority = Priority.High
        };
        await _alertRepository.AddAsync(alert);
        await _notificationService.NotifyAsync(alert);
    }
}
```

**✅ GOOD - Extracted Common Logic:**
```csharp
public async Task CheckTemperatureSensor(SensorId sensorId)
{
    var reading = await GetLatestReadingAsync(sensorId);
    if (reading.Value > TemperatureSafeLimit)
    {
        await RaiseCriticalAlertAsync(sensorId, "Temperature exceeds safe limit");
    }
}

public async Task CheckHumiditySensor(SensorId sensorId)
{
    var reading = await GetLatestReadingAsync(sensorId);
    if (reading.Value > HumiditySafeLimit)
    {
        await RaiseCriticalAlertAsync(sensorId, "Humidity exceeds safe limit");
    }
}

private async Task RaiseCriticalAlertAsync(SensorId sensorId, string message)
{
    var alert = CreateAlert(sensorId, AlertType.Critical, message, Priority.High);
    await SaveAndNotifyAlertAsync(alert);
}

private Alert CreateAlert(SensorId sensorId, AlertType type, string message, Priority priority)
{
    return new Alert
    {
        SensorId = sensorId,
        Type = type,
        Message = message,
        Timestamp = DateTime.UtcNow,
        Priority = priority
    };
}

private async Task SaveAndNotifyAlertAsync(Alert alert)
{
    await _alertRepository.AddAsync(alert);
    await _notificationService.NotifyAsync(alert);
}
```

---

## Summary

| Principle | Description |
|-----------|-------------|
| **Small** | 5-20 lines ideal, 30 max |
| **Do One Thing** | Single responsibility |
| **One Abstraction Level** | Don't mix high and low level |
| **Few Arguments** | 0-2 ideal, 3 max |
| **No Flag Arguments** | Split into separate methods |
| **Command-Query Separation** | Do or answer, not both |
| **No Side Effects** | Function does what name says |
| **Use Exceptions** | Not error codes |
| **DRY** | Don't repeat yourself |

**Remember:** Functions are the verbs of your program. Make them clear, focused, and expressive!
