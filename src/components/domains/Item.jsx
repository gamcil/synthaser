import { promiseOptions } from './CDD'
import AsyncSelect from 'react-select/async'


export const DomainItem = props => {

  // Mock up event since Select components do not
  // use event.target.value
  const onChange = domains => props.handleChange({
    target: {
      name: "domains",
      value: domains || []
    }
  })

  const defaultValues = domains => {
    if (!(domains && domains.length > 0))
      return null
    return domains.map(d => ({
      label: `${d.name} [${d.accession}]`,
      value: d.accession
    }))
  }

  return (
    <li>
      <button
        type="button"
        onClick={props.handleRemove}
      >
        Delete
      </button>
      <div className="rule-field">
        <label htmlFor="domain-type">Type:</label>
        <input
          id="domain-type"
          type="text"
          name="name"
          placeholder="Domain type..."
          value={props.data.name}
          onChange={props.handleChange}
          required
        />
      </div>
      <div className="rule-field">
        <label htmlFor="domain-families">Families:</label>
        <AsyncSelect
          id="domain-families"
          className="select"
          isMulti
          defaultValue={defaultValues(props.data.domains)}
          loadOptions={promiseOptions}
          onChange={onChange}
        />
      </div>
    </li>
  )
}
